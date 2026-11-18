"""Ensembl REST API client (GRCh38) for gene models, isoforms, CDS, Pfam domains, and paralogs."""

from __future__ import annotations

import urllib.error
import urllib.parse
from typing import Any, Optional

from crispr_preflight.cache import CachedHttpClient, NetworkFetchError
from crispr_preflight.models import ProteinDomain

ENSEMBL_BASE = "https://rest.ensembl.org"
EBI_PROTEINS_COORDS_BASE = "https://www.ebi.ac.uk/proteins/api/coordinates"
EBI_INTERPRO_PFAM_BASE = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot"
NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Primary GRCh38 (GCF_000001405.40) RefSeq chromosome accessions — 1-based inclusive coordinates identical to Ensembl GRCh38
GRCH38_NC_MAP: dict[str, str] = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X": "NC_000023.11",
    "Y": "NC_000024.10",
    "MT": "NC_012920.1",
}

CURATED_DOMAIN_DBS = {"Pfam", "UniProt", "Smart", "SuperFamily", "CDD", "PROSITE profiles"}
DOMAIN_DB_PRIORITY = {
    "Pfam": 0,
    "UniProt": 1,
    "Smart": 1,
    "CDD": 2,
    "PROSITE profiles": 3,
    "SuperFamily": 4,
}


def deduplicate_protein_domains(domains: list[ProteinDomain]) -> list[ProteinDomain]:
    """Collapse redundant cross-database annotations covering the same physical region.

    Ensembl protein_feature returns overlapping entries for the same structural fold
    across Pfam, Smart, CDD, PROSITE, and SuperFamily. We prioritize curated family
    models (Pfam > Smart > CDD > PROSITE > SuperFamily) and suppress lower-priority
    annotations whose amino-acid span overlaps >=50% of an already-retained domain,
    returning unique physical domains ordered N-terminus -> C-terminus.
    """
    ordered = sorted(
        domains,
        key=lambda d: (
            DOMAIN_DB_PRIORITY.get(d.source_db, 9),
            -(d.aa_end - d.aa_start),
            d.aa_start,
            d.accession,
        ),
    )
    kept: list[ProteinDomain] = []
    for cand in ordered:
        cand_len = max(1, cand.aa_end - cand.aa_start + 1)
        is_redundant = False
        for existing in kept:
            ov = min(cand.aa_end, existing.aa_end) - max(cand.aa_start, existing.aa_start) + 1
            if ov > 0:
                min_len = min(cand_len, max(1, existing.aa_end - existing.aa_start + 1))
                if (ov / min_len) >= 0.50:
                    is_redundant = True
                    break
        if not is_redundant:
            kept.append(cand)

    kept.sort(key=lambda d: (d.aa_start, d.aa_end, DOMAIN_DB_PRIORITY.get(d.source_db, 9)))
    return kept


class EnsemblClient:
    """Client for Ensembl REST endpoints with synonym fallback and batching."""

    def __init__(self, http: CachedHttpClient) -> None:
        self.http = http

    def fetch_ebi_coordinates(self, swissprot_id: str) -> Optional[dict[str, Any]]:
        """Fetch curated GRCh38 isoform genomic coordinates & protein features from EMBL-EBI Proteins API (~600 ms)."""
        if not swissprot_id:
            return None
        clean_acc = swissprot_id.strip().split("-")[0]
        url = f"{EBI_PROTEINS_COORDS_BASE}/{urllib.parse.quote(clean_acc)}"
        try:
            data = self.http.fetch_json(url, max_retries=1, timeout_override=6.0)
        except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
            return None
        if isinstance(data, dict) and data.get("gnCoordinate"):
            return data
        return None

    def fetch_interpro_pfam_domains(
        self,
        swissprot_id: str,
        ebi_coord: Optional[dict[str, Any]] = None,
    ) -> list[ProteinDomain]:
        """Fetch curated Pfam structural domains from EMBL-EBI InterPro API (~1.0s) + UniProt domain features."""
        if not swissprot_id:
            return []
        clean_acc = swissprot_id.strip().split("-")[0]
        url = f"{EBI_INTERPRO_PFAM_BASE}/{urllib.parse.quote(clean_acc)}/?page_size=50"
        domains: list[ProteinDomain] = []
        seen: set[tuple[str, str, int, int]] = set()

        try:
            ip_data = self.http.fetch_json(url, max_retries=1, timeout_override=6.0)
        except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
            ip_data = None

        if isinstance(ip_data, dict):
            for item in ip_data.get("results") or []:
                if not isinstance(item, dict):
                    continue
                meta = item.get("metadata") or {}
                acc = str(meta.get("accession") or "").strip()
                raw_name = meta.get("name")
                if isinstance(raw_name, dict):
                    desc = str(raw_name.get("name") or raw_name.get("short") or acc).strip()
                else:
                    desc = str(raw_name or acc).strip()
                for prot in item.get("proteins") or []:
                    if not isinstance(prot, dict):
                        continue
                    for loc in prot.get("entry_protein_locations") or []:
                        if not isinstance(loc, dict):
                            continue
                        for frag in loc.get("fragments") or []:
                            if not isinstance(frag, dict):
                                continue
                            start = int(frag.get("start") or 0)
                            end = int(frag.get("end") or 0)
                            if start <= 0 or end < start:
                                continue
                            dkey = ("Pfam", acc, start, end)
                            if dkey not in seen:
                                seen.add(dkey)
                                domains.append(
                                    ProteinDomain(
                                        source_db="Pfam",
                                        accession=acc or "Pfam",
                                        description=desc,
                                        aa_start=start,
                                        aa_end=end,
                                    )
                                )

        if isinstance(ebi_coord, dict):
            domains = self.merge_uniprot_domain_features(domains, ebi_coord)

        return deduplicate_protein_domains(domains)

    @staticmethod
    def merge_uniprot_domain_features(
        base_domains: list[ProteinDomain],
        ebi_coord: Optional[dict[str, Any]],
    ) -> list[ProteinDomain]:
        """Merge UniProt domain annotations from an EBI Proteins payload into an existing domain list."""
        if not isinstance(ebi_coord, dict):
            return deduplicate_protein_domains(base_domains)
        merged = list(base_domains)
        seen = {(d.source_db, d.description, d.aa_start, d.aa_end) for d in merged}
        for feat in ebi_coord.get("feature") or []:
            if not isinstance(feat, dict):
                continue
            if str(feat.get("type") or "").strip().lower() != "domain":
                continue
            desc = str(feat.get("description") or "Conserved domain").strip()
            loc = feat.get("location") or {}
            b_pos = ((loc.get("begin") or {}).get("position")) if isinstance(loc, dict) else None
            e_pos = ((loc.get("end") or {}).get("position")) if isinstance(loc, dict) else None
            if b_pos is None or e_pos is None:
                continue
            start, end = min(int(b_pos), int(e_pos)), max(int(b_pos), int(e_pos))
            if start <= 0 or end < start:
                continue
            dkey = ("UniProt", desc, start, end)
            if dkey not in seen:
                seen.add(dkey)
                merged.append(
                    ProteinDomain(
                        source_db="UniProt",
                        accession="UniProtKB",
                        description=desc,
                        aa_start=start,
                        aa_end=end,
                    )
                )
        return deduplicate_protein_domains(merged)

    @staticmethod
    def _build_gene_from_ebi_coords(
        ensembl_id: str,
        gene_symbol: str,
        approved_name: str,
        gc_list: list[dict[str, Any]],
        canon_seq_len: int,
    ) -> Optional[dict[str, Any]]:
        """Synthesize an Ensembl-compatible GRCh38 Gene dict directly from EBI Proteins gnCoordinate."""
        txs: list[dict[str, Any]] = []
        canon_assigned = False
        chrom_out = ""
        strand_out = 1

        for idx, gc in enumerate(gc_list, start=1):
            enst = str(gc.get("ensemblTranscriptId") or "").split(".")[0]
            ensp = str(gc.get("ensemblTranslationId") or "").split(".")[0]
            gl = gc.get("genomicLocation") or {}
            if not enst or not isinstance(gl, dict):
                continue
            chrom = str(gl.get("chromosome") or "")
            if not chrom or chrom.startswith("CHR_"):
                continue
            strand = -1 if gl.get("reverseStrand") else 1
            raw_exons = gl.get("exon") or []
            exons_out: list[dict[str, Any]] = []
            max_prot_aa = 0
            for ex_idx, ex in enumerate(raw_exons, start=1):
                if not isinstance(ex, dict):
                    continue
                gloc = ex.get("genomeLocation") or {}
                b = (gloc.get("begin") or {}).get("position")
                e = (gloc.get("end") or {}).get("position")
                if b is None or e is None:
                    continue
                s_pos, e_pos = min(int(b), int(e)), max(int(b), int(e))
                exons_out.append(
                    {
                        "id": str(ex.get("id") or f"{enst}.E{ex_idx}"),
                        "start": s_pos,
                        "end": e_pos,
                        "strand": strand,
                    }
                )
                ploc = ex.get("proteinLocation") or {}
                pe = (ploc.get("end") or {}).get("position")
                if pe is not None:
                    max_prot_aa = max(max_prot_aa, int(pe))

            if not exons_out:
                continue

            tx_start = min(e["start"] for e in exons_out)
            tx_end = max(e["end"] for e in exons_out)
            cds_bp = sum(e["end"] - e["start"] + 1 for e in exons_out)
            prot_aa = max_prot_aa if max_prot_aa > 0 else max(1, cds_bp // 3)
            is_canon = False
            if not canon_assigned and canon_seq_len > 0 and prot_aa == canon_seq_len:
                is_canon = True
                canon_assigned = True

            if not chrom_out:
                chrom_out = chrom
                strand_out = strand

            txs.append(
                {
                    "id": enst,
                    "object_type": "Transcript",
                    "display_name": f"{gene_symbol}-20{idx}",
                    "biotype": "protein_coding",
                    "is_canonical": 1 if is_canon else 0,
                    "is_mane_select": is_canon,
                    "gencode_primary": 1,
                    "seq_region_name": chrom,
                    "strand": strand,
                    "start": tx_start,
                    "end": tx_end,
                    "Translation": {
                        "id": ensp or f"{enst}_P",
                        "start": tx_start,
                        "end": tx_end,
                        "length": prot_aa,
                    },
                    "Exon": exons_out,
                }
            )

        if not txs:
            return None
        if not canon_assigned:
            txs[0]["is_canonical"] = 1
            txs[0]["is_mane_select"] = True

        return {
            "id": ensembl_id,
            "object_type": "Gene",
            "display_name": gene_symbol,
            "description": approved_name,
            "seq_region_name": chrom_out,
            "strand": strand_out,
            "start": min(int(t["start"]) for t in txs),
            "end": max(int(t["end"]) for t in txs),
            "Transcript": txs,
        }

    def lookup_gene_fast(
        self,
        ensembl_id: str,
        gene_symbol: str,
        approved_name: str = "",
        ebi_coord: Optional[dict[str, Any]] = None,
    ) -> tuple[dict[str, Any], bool]:
        """Resolve GRCh38 gene & curated isoforms in ~1s instead of expanding 80+ computational transcripts (~37s)."""
        if isinstance(ebi_coord, dict) and ebi_coord.get("gnCoordinate"):
            canon_seq_len = len(str(ebi_coord.get("sequence") or ""))
            gc_matches: list[dict[str, Any]] = []
            for gc in ebi_coord.get("gnCoordinate") or []:
                if not isinstance(gc, dict):
                    continue
                gc_gene = str(gc.get("ensemblGeneId") or "").split(".")[0]
                gl = gc.get("genomicLocation") or {}
                asm = str(gl.get("assemblyName") or "GRCh38").upper()
                chrom = str(gl.get("chromosome") or "")
                if gc_gene == ensembl_id and "GRCH38" in asm and chrom and not chrom.startswith("CHR_"):
                    gc_matches.append(gc)

            if gc_matches:
                # Build GRCh38 gene & Swiss-Prot curated isoform models directly in <0.1 ms
                # without waiting on rest.ensembl.org's overloaded /lookup/ endpoint.
                synth = self._build_gene_from_ebi_coords(
                    ensembl_id=ensembl_id,
                    gene_symbol=gene_symbol,
                    approved_name=approved_name,
                    gc_list=gc_matches,
                    canon_seq_len=canon_seq_len,
                )
                if synth is not None:
                    return synth, True

        return self.lookup_gene_expanded(ensembl_id), False

    def lookup_gene_expanded(self, gene_query: str) -> dict[str, Any]:
        """Resolve a gene symbol, alias (e.g. BRG1 -> SMARCA4), or ENSG ID with expanded transcripts."""
        cleaned = gene_query.strip()
        if cleaned.upper().startswith("ENSG"):
            base_id = cleaned.upper().split(".")[0]
            url = f"{ENSEMBL_BASE}/lookup/id/{urllib.parse.quote(base_id)}?expand=1;mane=1"
            data = self.http.fetch_json(url)
            if isinstance(data, dict) and data.get("object_type") == "Gene":
                return data
            raise ValueError(f"Ensembl ID '{cleaned}' was not found as a human gene.")

        # 1. Fast lightweight symbol -> ENSG resolution (expand=0 avoids Ensembl symbol+expand timeout)
        ensg_id: Optional[str] = None
        sym_url = f"{ENSEMBL_BASE}/lookup/symbol/human/{urllib.parse.quote(cleaned)}?expand=0"
        try:
            sym_data = self.http.fetch_json(sym_url)
            if isinstance(sym_data, dict) and str(sym_data.get("id", "")).startswith("ENSG"):
                ensg_id = str(sym_data["id"])
        except urllib.error.HTTPError as exc:
            if exc.code not in (400, 404):
                raise

        if ensg_id:
            id_url = f"{ENSEMBL_BASE}/lookup/id/{ensg_id}?expand=1;mane=1"
            data = self.http.fetch_json(id_url)
            if isinstance(data, dict) and data.get("object_type") == "Gene":
                return data

        # 2. Synonym / alias fallback (e.g. BRG1, ERK2, HER2, P53)
        xref_url = f"{ENSEMBL_BASE}/xrefs/symbol/human/{urllib.parse.quote(cleaned)}?object_type=gene"
        xrefs = self.http.fetch_json(xref_url) or []
        ensg_candidates = [
            x["id"]
            for x in xrefs
            if isinstance(x, dict) and str(x.get("id", "")).startswith("ENSG")
        ]
        for cand_id in ensg_candidates:
            id_url = f"{ENSEMBL_BASE}/lookup/id/{cand_id}?expand=1;mane=1"
            try:
                gene_obj = self.http.fetch_json(id_url)
                if isinstance(gene_obj, dict) and not str(gene_obj.get("seq_region_name", "")).startswith("CHR_"):
                    return gene_obj
            except (urllib.error.HTTPError, NetworkFetchError):
                continue

        if ensg_candidates:
            id_url = f"{ENSEMBL_BASE}/lookup/id/{ensg_candidates[0]}?expand=1;mane=1"
            return self.http.fetch_json(id_url)

        raise ValueError(
            f"Could not resolve human gene symbol or alias '{gene_query}' in Ensembl GRCh38."
        )

    def fetch_cds_sequence(self, transcript_id: str) -> str:
        """Fetch stitched 5'->3' coding nucleotide sequence (ATG...STOP) for an ENST transcript."""
        url = f"{ENSEMBL_BASE}/sequence/id/{urllib.parse.quote(transcript_id)}?type=cds"
        raw = self.http.fetch_bytes(url, headers={"Accept": "text/plain"})
        return raw.decode("utf-8").strip().upper()

    def fetch_protein_domains(self, protein_id: str) -> list[ProteinDomain]:
        """Fetch conserved structural/functional domains mapped onto ENSP amino-acid coordinates."""
        if not protein_id:
            return []
        url = f"{ENSEMBL_BASE}/overlap/translation/{urllib.parse.quote(protein_id)}?feature=protein_feature"
        try:
            items = self.http.fetch_json(url) or []
        except (urllib.error.HTTPError, NetworkFetchError):
            return []

        domains: list[ProteinDomain] = []
        seen: set[tuple[str, str, int, int]] = set()
        filtered = [
            it for it in items
            if isinstance(it, dict) and it.get("type") in CURATED_DOMAIN_DBS
        ]

        for item in filtered:
            db = str(item.get("type", ""))
            acc = str(item.get("id", ""))
            desc = str(item.get("description") or acc).strip()
            start = int(item.get("start", 0))
            end = int(item.get("end", 0))
            if start <= 0 or end < start:
                continue
            dedup_key = (db, acc, start, end)
            if dedup_key in seen:
                continue
            seen.add(dedup_key)
            domains.append(
                ProteinDomain(
                    source_db=db,
                    accession=acc,
                    description=desc,
                    aa_start=start,
                    aa_end=end,
                )
            )
        return deduplicate_protein_domains(domains)

    def fetch_human_paralogs(self, ensembl_id: str, min_identity: float = 30.0) -> list[dict[str, Any]]:
        """Fetch human within-species paralogs from Ensembl Compara and resolve their gene symbols."""
        url = (
            f"{ENSEMBL_BASE}/homology/id/human/"
            f"{urllib.parse.quote(ensembl_id)}?type=paralogues;sequence=none"
        )
        try:
            resp = self.http.fetch_json(url) or {}
        except (urllib.error.HTTPError, NetworkFetchError):
            return []

        data_list = resp.get("data", [])
        if not data_list:
            return []

        homologies = data_list[0].get("homologies", [])
        raw_paralogs: list[dict[str, Any]] = []
        target_ids: list[str] = []

        for hom in homologies:
            target = hom.get("target", {})
            source = hom.get("source", {})
            if target.get("species") != "homo_sapiens":
                continue
            t_id = str(target.get("id", ""))
            if not t_id or t_id == ensembl_id:
                continue
            q_ident = float(source.get("perc_id") or 0.0)
            t_ident = float(target.get("perc_id") or 0.0)
            if max(q_ident, t_ident) < min_identity:
                continue
            raw_paralogs.append(
                {
                    "ensembl_id": t_id,
                    "percent_identity": round(q_ident, 1),
                    "target_percent_identity": round(t_ident, 1),
                    "ancestor_taxonomy": hom.get("taxonomy_level"),
                }
            )
            target_ids.append(t_id)

        if not raw_paralogs:
            return []

        # Batch-resolve ENSG -> display_name in a single POST call
        symbol_map = self._batch_lookup_symbols(target_ids[:40])
        for entry in raw_paralogs:
            entry["gene_symbol"] = symbol_map.get(entry["ensembl_id"], entry["ensembl_id"])

        # Deduplicate by symbol and sort by query sequence identity descending
        deduped: dict[str, dict[str, Any]] = {}
        for entry in sorted(raw_paralogs, key=lambda x: x["percent_identity"], reverse=True):
            sym = entry["gene_symbol"]
            if sym not in deduped:
                deduped[sym] = entry
        return list(deduped.values())

    def _batch_lookup_symbols(self, ensg_ids: list[str]) -> dict[str, str]:
        """Resolve up to 50 Ensembl IDs to HGNC display symbols in one POST call."""
        if not ensg_ids:
            return {}
        url = f"{ENSEMBL_BASE}/lookup/id?expand=0"
        try:
            resp = self.http.fetch_json(url, method="POST", json_body={"ids": ensg_ids}) or {}
        except (urllib.error.HTTPError, NetworkFetchError):
            return {}
        out: dict[str, str] = {}
        for k, val in resp.items():
            if isinstance(val, dict) and val.get("display_name"):
                out[k] = str(val["display_name"])
        return out

    def _fetch_flanked_via_ncbi(
        self,
        nc_acc: str,
        intervals: list[tuple[int, int]],
        pad: int,
        max_cluster_span: int = 250_000,
        max_clusters: int = 3,
    ) -> Optional[dict[tuple[int, int], Optional[str]]]:
        """Slice exact 1-based GRCh38 forward-strand windows from NCBI E-utilities (~550 ms vs 7-12s on Ensembl)."""
        ordered = sorted(set(intervals), key=lambda iv: (max(1, iv[0] - pad), iv[1] + pad))
        clusters: list[tuple[int, int, list[tuple[int, int]]]] = []
        for s, e in ordered:
            ps, pe = max(1, s - pad), e + pad
            if (pe - ps + 1) > max_cluster_span:
                return None
            if not clusters:
                clusters.append((ps, pe, [(s, e)]))
            else:
                c_s, c_e, c_list = clusters[-1]
                new_e = max(c_e, pe)
                if (new_e - c_s + 1) <= max_cluster_span:
                    clusters[-1] = (c_s, new_e, c_list + [(s, e)])
                else:
                    clusters.append((ps, pe, [(s, e)]))
                    if len(clusters) > max_clusters:
                        return None

        out: dict[tuple[int, int], Optional[str]] = {}
        for c_s, c_e, c_list in clusters:
            url = (
                f"{NCBI_EFETCH_URL}?db=nuccore&id={urllib.parse.quote(nc_acc)}"
                f"&seq_start={c_s}&seq_stop={c_e}&strand=1&rettype=fasta&retmode=text"
            )
            try:
                raw = self.http.fetch_bytes(
                    url,
                    headers={"Accept": "text/plain"},
                    max_retries=1,
                    timeout_override=6.5,
                )
            except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
                return None
            seq = "".join(
                line.strip()
                for line in raw.decode("utf-8", errors="replace").splitlines()
                if line and not line.startswith(">")
            ).upper()
            expected_len = c_e - c_s + 1
            if len(seq) != expected_len:
                return None
            for s, e in c_list:
                ps, pe = max(1, s - pad), e + pad
                sub = seq[ps - c_s : pe - c_s + 1]
                if len(sub) != (pe - ps + 1):
                    return None
                out[(s, e)] = sub
        return out

    def fetch_flanked_regions(
        self,
        chrom: str,
        intervals: list[tuple[int, int]],
        pad: int = 25,
        nc_chrom_map: Optional[dict[str, str]] = None,
        ensembl_species: str = "human",
    ) -> dict[tuple[int, int], Optional[str]]:
        """Batch-fetch +1 forward-strand genomic sequences for [start-pad, end+pad] windows."""
        if not intervals:
            return {}

        active_nc_map = nc_chrom_map if nc_chrom_map is not None else GRCH38_NC_MAP
        clean_chr = str(chrom).replace("chr", "").upper()
        nc_acc = active_nc_map.get(clean_chr) or active_nc_map.get(str(chrom).replace("chr", ""))
        if nc_acc:
            fast_res = self._fetch_flanked_via_ncbi(nc_acc, intervals, pad)
            if fast_res is not None and all(fast_res.get(iv) for iv in intervals):
                return fast_res

        result: dict[tuple[int, int], Optional[str]] = {}
        # Fallback: Ensembl POST /sequence/region accepts up to 50 regions per batch
        batch_size = 45
        for idx in range(0, len(intervals), batch_size):
            chunk = intervals[idx : idx + batch_size]
            region_specs: list[str] = []
            spec_to_interval: dict[int, tuple[int, int]] = {}
            for j, (s, e) in enumerate(chunk):
                p_start = max(1, s - pad)
                p_end = e + pad
                region_specs.append(f"{chrom}:{p_start}..{p_end}:1")
                spec_to_interval[j] = (s, e)

            url = f"{ENSEMBL_BASE}/sequence/region/{urllib.parse.quote(ensembl_species)}"
            try:
                resp = self.http.fetch_json(
                    url,
                    method="POST",
                    json_body={"regions": region_specs},
                )
                if isinstance(resp, list):
                    for j, item in enumerate(resp):
                        orig = spec_to_interval.get(j)
                        if orig is not None and isinstance(item, dict):
                            result[orig] = str(item.get("seq", "")).upper() or None
            except (urllib.error.HTTPError, NetworkFetchError):
                for s, e in chunk:
                    result[(s, e)] = None
        return result
