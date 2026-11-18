"""Fast clients for MyGene.info TSS coordinates, cBioPortal CCLE cell-line genetics, and MyVariant hg38 SNPs."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import threading
import urllib.error
import urllib.parse
from typing import Any, Optional

from crispr_preflight.cache import CachedHttpClient, NetworkFetchError
from crispr_preflight.models import CellLineAlterationProfile

MYGENE_GENE_URL = "https://mygene.info/v3/gene"
MYVARIANT_QUERY_URL = "https://myvariant.info/v1/query"
CBIOPORTAL_BASE = "https://www.cbioportal.org/api"

TRUNCATING_MUTATION_TYPES = {
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Nonsense_Mutation",
    "Splice_Site",
    "Nonstop_Mutation",
    "Translation_Start_Site",
}

CNA_LABELS: dict[int, str] = {
    -2: "Homozygous Deletion (GISTIC -2)",
    -1: "Hemizygous Loss (GISTIC -1)",
    0: "Copy-Neutral (GISTIC 0)",
    1: "Low-Level Copy Gain (GISTIC +1)",
    2: "High-Level Amplification (GISTIC +2)",
}

_CCLE_LOCK = threading.Lock()
_CCLE_SAMPLES_CACHE: Optional[dict[str, str]] = None


def _norm_cell_token(name: str) -> str:
    return "".join(ch for ch in name.strip().upper() if ch.isalnum())


class CellLineAndVariantClient:
    """Queries MyGene.info (TSS), cBioPortal CCLE (cell-line LoF/CNA), and MyVariant.info (common SNPs)."""

    def __init__(self, http: CachedHttpClient) -> None:
        self.http = http

    def fetch_tss_from_mygene(
        self,
        ensembl_id: str,
        chrom: str = "",
        strand: Optional[int] = None,
    ) -> Optional[int]:
        """Resolve the 1-based GRCh38 5' Transcription Start Site (TSS) of the primary RefSeq transcript (~700 ms)."""
        if not ensembl_id:
            return None
        url = f"{MYGENE_GENE_URL}/{urllib.parse.quote(ensembl_id)}?fields=exons,genomic_pos"
        try:
            data = self.http.fetch_json(url, max_retries=0, timeout_override=3.5)
        except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
            return None
        if not isinstance(data, dict):
            return None

        clean_chr = chrom.replace("chr", "").upper()
        raw_exons = data.get("exons")
        tx_list: list[dict[str, Any]] = []
        if isinstance(raw_exons, list):
            tx_list = [x for x in raw_exons if isinstance(x, dict)]
        elif isinstance(raw_exons, dict):
            tx_list = [v for v in raw_exons.values() if isinstance(v, dict)]

        valid_txs: list[tuple[int, int, dict[str, Any]]] = []
        for tx in tx_list:
            tx_chr = str(tx.get("chr") or "").replace("chr", "").upper()
            if clean_chr and tx_chr and tx_chr != clean_chr:
                continue
            if "_" in tx_chr or tx_chr.startswith("CHR_"):
                continue
            tx_strand = int(tx.get("strand") or 1)
            if strand is not None and tx_strand != strand:
                continue
            tx_s = tx.get("txstart")
            tx_e = tx.get("txend")
            if tx_s is None or tx_e is None:
                continue
            acc = str(tx.get("transcript") or "")
            is_nm = 0 if acc.startswith("NM_") else 1
            cds_span = abs(int(tx.get("cdsend") or 0) - int(tx.get("cdsstart") or 0))
            valid_txs.append((is_nm, -cds_span, tx))

        if valid_txs:
            valid_txs.sort(key=lambda item: (item[0], item[1], str(item[2].get("transcript") or "")))
            best = valid_txs[0][2]
            b_strand = int(best.get("strand") or 1)
            # RefSeq/UCSC txstart is 0-based half-open [txstart, txend)
            return (int(best["txstart"]) + 1) if b_strand >= 0 else int(best["txend"])

        gp = data.get("genomic_pos")
        if isinstance(gp, list):
            gp = next(
                (g for g in gp if isinstance(g, dict) and not str(g.get("chr", "")).startswith("CHR_")),
                None,
            )
        if isinstance(gp, dict) and gp.get("start") is not None and gp.get("end") is not None:
            g_strand = int(gp.get("strand") or strand or 1)
            return int(gp["start"]) if g_strand >= 0 else int(gp["end"])
        return None

    def _get_ccle_sample_map(self) -> dict[str, str]:
        """Return {NORMALIZED_CELL_NAME: 'SAMPLE_ID'} across 1,739 CCLE cell lines (cached 7d on disk + RAM)."""
        global _CCLE_SAMPLES_CACHE
        with _CCLE_LOCK:
            if _CCLE_SAMPLES_CACHE is not None:
                return _CCLE_SAMPLES_CACHE
        url = f"{CBIOPORTAL_BASE}/studies/ccle_broad_2019/samples?projection=SUMMARY&pageSize=2500"
        try:
            items = self.http.fetch_json(url, max_retries=0, timeout_override=4.0) or []
        except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
            return {}
        mapping: dict[str, str] = {}
        if isinstance(items, list):
            for it in items:
                if not isinstance(it, dict) or not it.get("sampleId"):
                    continue
                sid = str(it["sampleId"])
                prefix = sid.split("_", 1)[0]
                norm = _norm_cell_token(prefix)
                if norm and norm not in mapping:
                    mapping[norm] = sid
        if mapping:
            with _CCLE_LOCK:
                _CCLE_SAMPLES_CACHE = mapping
        return mapping

    def resolve_ccle_sample_id(self, cell_line_query: str) -> Optional[str]:
        norm = _norm_cell_token(cell_line_query)
        if not norm:
            return None
        smap = self._get_ccle_sample_map()
        if not smap:
            return None
        if norm in smap:
            return smap[norm]
        # Strip trailing T/S subclone suffix (e.g. HEK293T -> HEK293) or check prefix match
        for k, sid in smap.items():
            if len(norm) >= 4 and (k == norm or sid.upper().startswith(norm + "_")):
                return sid
        return None

    def fetch_ccle_alterations(
        self,
        gene_symbol: str,
        cell_line_query: Optional[str],
    ) -> Optional[CellLineAlterationProfile]:
        """Fetch endogenous somatic mutations and GISTIC copy-number status in the requested CCLE cell line."""
        if not cell_line_query or not cell_line_query.strip():
            return None
        profile = CellLineAlterationProfile(requested_cell_line=cell_line_query.strip())
        sample_id = self.resolve_ccle_sample_id(cell_line_query)
        if not sample_id:
            return profile
        profile.matched_ccle_sample = sample_id

        gene_url = f"{CBIOPORTAL_BASE}/genes/{urllib.parse.quote(gene_symbol.strip().upper())}"
        try:
            g_obj = self.http.fetch_json(gene_url, max_retries=0, timeout_override=3.2)
        except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
            return profile
        if not isinstance(g_obj, dict) or not g_obj.get("entrezGeneId"):
            return profile
        entrez_id = int(g_obj["entrezGeneId"])

        mut_url = f"{CBIOPORTAL_BASE}/molecular-profiles/ccle_broad_2019_mutations/mutations/fetch?projection=DETAILED"
        cna_url = f"{CBIOPORTAL_BASE}/molecular-profiles/ccle_broad_2019_cna/molecular-data/fetch"
        payload = {"entrezGeneIds": [entrez_id], "sampleIds": [sample_id]}

        def _fetch_muts() -> list[dict[str, Any]]:
            try:
                res = self.http.fetch_json(
                    mut_url, method="POST", json_body=payload, max_retries=0, timeout_override=3.5
                )
                return res if isinstance(res, list) else []
            except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
                return []

        def _fetch_cna() -> list[dict[str, Any]]:
            try:
                res = self.http.fetch_json(
                    cna_url, method="POST", json_body=payload, max_retries=0, timeout_override=3.5
                )
                return res if isinstance(res, list) else []
            except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
                return []

        with ThreadPoolExecutor(max_workers=2) as pool:
            f_m = pool.submit(_fetch_muts)
            f_c = pool.submit(_fetch_cna)
            mut_rows = f_m.result()
            cna_rows = f_c.result()

        seen_muts: set[str] = set()
        for m in mut_rows:
            if not isinstance(m, dict):
                continue
            mtype = str(m.get("mutationType") or "Variant")
            if mtype == "Silent":
                continue
            pchange = str(m.get("proteinChange") or "").strip()
            p_lbl = f"p.{pchange}" if (pchange and not pchange.startswith("p.")) else (pchange or "coding variant")
            vaf = m.get("variantAlleleFrequency")
            vaf_str = f", VAF {int(round(float(vaf) * 100))}%" if isinstance(vaf, (int, float)) and vaf > 0 else ""
            desc = f"{p_lbl} ({mtype}{vaf_str})"
            if desc not in seen_muts:
                seen_muts.add(desc)
                profile.mutations.append(desc)
            if mtype in TRUNCATING_MUTATION_TYPES:
                profile.has_truncating_lof = True

        if cna_rows and isinstance(cna_rows[0], dict) and cna_rows[0].get("value") is not None:
            try:
                cval = int(round(float(cna_rows[0]["value"])))
                profile.cna_gistic = cval
                profile.cna_label = CNA_LABELS.get(cval, f"GISTIC {cval}")
            except (ValueError, TypeError):
                pass

        return profile

    def check_guides_common_snps(
        self,
        guide_specs: list[dict[str, Any]],
        min_af: float = 0.01,
    ) -> dict[str, dict[str, Any]]:
        """Batch-check up to 20 guide+PAM 23-bp GRCh38 windows for common human SNPs (gnomAD AF >= 1%)."""
        out: dict[str, dict[str, Any]] = {}
        valid_specs: list[dict[str, Any]] = []
        clauses: list[str] = []

        for sp in guide_specs[:20]:
            gid = str(sp.get("guide_id") or "")
            chrom = str(sp.get("chrom") or "").replace("chr", "")
            s = sp.get("start")
            e = sp.get("end")
            strand = str(sp.get("strand") or "+")
            if not gid or not chrom or s is None or e is None:
                continue
            s_i, e_i = min(int(s), int(e)), max(int(s), int(e))
            pam_orient = str(sp.get("pam_orientation") or "3prime")
            pam_len = int(sp.get("pam_len") or 3)
            spacer_len = int(sp.get("spacer_len") or max(20, (e_i - s_i + 1) - pam_len))
            out[gid] = {"status": "CLEAN", "summary": "No common SNPs (MAF >=1%)", "variants": []}
            valid_specs.append(
                {
                    "guide_id": gid,
                    "chrom": chrom,
                    "start": s_i,
                    "end": e_i,
                    "strand": strand,
                    "pam_orientation": pam_orient,
                    "pam_len": pam_len,
                    "spacer_len": spacer_len,
                }
            )
            clause = f"chr{chrom}:{s_i}-{e_i}"
            if clause not in clauses:
                clauses.append(clause)

        if not clauses:
            return out

        q_str = (
            f"({' OR '.join(clauses)}) AND "
            f"(gnomad_genome.af.af:>{min_af} OR gnomad_exome.af.af:>{min_af})"
        )
        url = (
            f"{MYVARIANT_QUERY_URL}?q={urllib.parse.quote(q_str)}"
            "&assembly=hg38&fields=dbsnp.rsid,gnomad_genome.af.af,gnomad_exome.af.af,hg38&size=50"
        )
        try:
            resp = self.http.fetch_json(url, max_retries=0, timeout_override=5.0)
        except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
            for gid in out:
                out[gid] = {"status": "UNAVAILABLE", "summary": "SNP lookup timed out", "variants": []}
            return out

        hits = (resp or {}).get("hits") if isinstance(resp, dict) else []
        for h in hits or []:
            if not isinstance(h, dict):
                continue
            hg38 = h.get("hg38") or {}
            if isinstance(hg38, list) and hg38:
                hg38 = hg38[0]
            vs = hg38.get("start")
            ve = hg38.get("end")
            if vs is None:
                continue
            vs_i = int(vs)
            ve_i = int(ve) if ve is not None else vs_i
            vid_chr = str(h.get("_id") or "").split(":", 1)[0].replace("chr", "")

            rsid = str(((h.get("dbsnp") or {}).get("rsid")) or h.get("_id") or "variant")
            af_g = float((((h.get("gnomad_genome") or {}).get("af") or {}).get("af")) or 0.0)
            af_e = float((((h.get("gnomad_exome") or {}).get("af") or {}).get("af")) or 0.0)
            af = max(af_g, af_e)
            if af < min_af:
                continue
            af_pct = round(af * 100.0, 1)

            for sp in valid_specs:
                if vid_chr and sp["chrom"] != vid_chr:
                    continue
                if max(vs_i, sp["start"]) <= min(ve_i, sp["end"]):
                    # Classify exact zone along the 5'->3' target strand ([spacer][PAM] or [PAM][spacer])
                    win_max = max(0, sp["end"] - sp["start"])
                    pos_5to3 = (
                        max(0, min(win_max, vs_i - sp["start"]))
                        if sp["strand"] == "+"
                        else max(0, min(win_max, sp["end"] - vs_i))
                    )
                    if sp["pam_orientation"] == "5prime":
                        p_len = sp["pam_len"]
                        zone = "PAM" if pos_5to3 < p_len else ("Seed" if pos_5to3 < p_len + 8 else "Distal")
                    else:
                        s_len = sp["spacer_len"]
                        zone = "PAM" if pos_5to3 >= s_len else ("Seed" if pos_5to3 >= s_len - 10 else "Distal")

                    entry = out[sp["guide_id"]]
                    v_label = f"{rsid} (MAF {af_pct}% in {zone})"
                    if v_label not in entry["variants"]:
                        entry["variants"].append(v_label)
                    entry["status"] = "HIGH" if zone in ("PAM", "Seed") else "WARN"
                    entry["summary"] = "; ".join(entry["variants"])

        return out
