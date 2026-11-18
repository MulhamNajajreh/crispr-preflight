"""Fast multi-species gene model, RefSeq isoform, Alliance paralog, and human-ortholog resolver.

Supports Mouse (GRCm39), Rat (GRCr8), Zebrafish (GRCz11), Fruit Fly (BDGP6.46), and
C. elegans (WBcel235) in ~1.6-2.4s cold without relying on slow Ensembl REST expansions.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import urllib.error
import urllib.parse
from typing import Any, Optional

from crispr_preflight.cache import CachedHttpClient, NetworkFetchError
from crispr_preflight.clients.ensembl import EnsemblClient
from crispr_preflight.species import SpeciesConfig

MYGENE_QUERY_URL = "https://mygene.info/v3/query"
NCBI_DATASETS_GENE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id"
ALLIANCE_GENE_URL = "https://www.alliancegenome.org/api/gene"

_CANON_ENS_PREFIX: dict[str, str] = {
    "mouse": "ENSMUSG",
    "rat": "ENSRNOG",
    "zebrafish": "ENSDARG",
    "fly": "FBgn",
    "worm": "WBGene",
}


class MultiSpeciesClient:
    """Resolves non-human model organism gene models, exons, paralogs, and human orthologs."""

    def __init__(self, http: CachedHttpClient, ensembl: EnsemblClient) -> None:
        self.http = http
        self.ensembl = ensembl

    @staticmethod
    def _extract_ensembl_gene_id(raw_ens: Any, species_id: str) -> Optional[str]:
        pref = _CANON_ENS_PREFIX.get(species_id, "ENS")
        candidates: list[str] = []
        if isinstance(raw_ens, dict) and raw_ens.get("gene"):
            candidates.append(str(raw_ens["gene"]))
        elif isinstance(raw_ens, list):
            for item in raw_ens:
                if isinstance(item, dict) and item.get("gene"):
                    candidates.append(str(item["gene"]))
        for c in candidates:
            if c.startswith(pref):
                return c
        return candidates[0] if candidates else None

    @staticmethod
    def _extract_swissprot(raw_up: Any) -> Optional[str]:
        if not isinstance(raw_up, dict):
            return None
        sp = raw_up.get("Swiss-Prot")
        if isinstance(sp, list) and sp:
            return str(sp[0]).split("-")[0]
        if isinstance(sp, str) and sp.strip():
            return sp.strip().split("-")[0]
        return None

    @staticmethod
    def _extract_mod_curie(hit: dict[str, Any], species_id: str) -> Optional[str]:
        if species_id == "mouse" and hit.get("MGI"):
            val = str(hit["MGI"])
            return val if val.startswith("MGI:") else f"MGI:{val}"
        if species_id == "rat" and hit.get("RGD"):
            val = str(hit["RGD"])
            return val if val.startswith("RGD:") else f"RGD:{val}"
        if species_id == "zebrafish" and hit.get("ZFIN"):
            val = str(hit["ZFIN"])
            return val if val.startswith("ZFIN:") else f"ZFIN:{val}"
        if species_id == "fly" and hit.get("FLYBASE"):
            val = str(hit["FLYBASE"])
            return val if val.startswith("FB:") else f"FB:{val}"
        if species_id == "worm" and hit.get("WormBase"):
            val = str(hit["WormBase"])
            return val if val.startswith("WB:") else f"WB:{val}"
        return None

    def query_mygene_hit(self, gene_query: str, species_cfg: SpeciesConfig) -> dict[str, Any]:
        """Resolve a symbol, alias, or gene ID in the target organism via MyGene.info (~500-700 ms)."""
        cleaned = gene_query.strip()
        fields = "symbol,name,alias,entrezgene,ensembl.gene,uniprot,MGI,RGD,ZFIN,FLYBASE,WormBase,homologene,exons,genomic_pos"
        url = (
            f"{MYGENE_QUERY_URL}?q={urllib.parse.quote(cleaned)}"
            f"&species={species_cfg.taxid}&fields={fields}&size=10"
        )
        data = self.http.fetch_json(url, max_retries=1, timeout_override=6.0)
        hits = (data or {}).get("hits") if isinstance(data, dict) else []
        if not hits:
            raise ValueError(
                f"Could not resolve '{gene_query}' in {species_cfg.common_name} ({species_cfg.scientific_name}, {species_cfg.assembly})."
            )

        q_low = cleaned.lower()
        # Rank 0: exact symbol match; Rank 1: exact Ensembl/MOD ID match; Rank 2: exact alias match; Rank 3: top MyGene hit
        best_hit = hits[0]
        best_rank = 9
        for h in hits:
            if not isinstance(h, dict):
                continue
            sym = str(h.get("symbol") or "")
            ens_id = self._extract_ensembl_gene_id(h.get("ensembl"), species_cfg.species_id) or ""
            mod_curie = self._extract_mod_curie(h, species_cfg.species_id) or ""
            aliases_raw = h.get("alias")
            aliases = (
                [str(a).lower() for a in aliases_raw]
                if isinstance(aliases_raw, list)
                else ([str(aliases_raw).lower()] if aliases_raw else [])
            )
            if sym.lower() == q_low:
                rank = 0
            elif q_low in (ens_id.lower(), mod_curie.lower(), mod_curie.split(":")[-1].lower()):
                rank = 1
            elif q_low in aliases:
                rank = 2
            else:
                rank = 3
            # Prefer hits that carry coding exons or Entrez IDs
            has_exons = 0 if (h.get("exons") or h.get("entrezgene")) else 1
            if (rank, has_exons) < (best_rank, 0 if (best_hit.get("exons") or best_hit.get("entrezgene")) else 1):
                best_rank = rank
                best_hit = h
                if rank == 0 and has_exons == 0:
                    break

        return best_hit

    @staticmethod
    def build_gene_from_mygene_exons(
        hit: dict[str, Any],
        gene_id: str,
        gene_symbol: str,
        description: str,
        species_cfg: SpeciesConfig,
    ) -> tuple[Optional[dict[str, Any]], Optional[int]]:
        """Build an Ensembl-compatible Gene dict + 1-based 5' TSS from MyGene RefSeq `exons`."""
        raw_exons = hit.get("exons")
        tx_items: list[dict[str, Any]] = []
        if isinstance(raw_exons, list):
            tx_items = [x for x in raw_exons if isinstance(x, dict)]
        elif isinstance(raw_exons, dict):
            tx_items = [v for v in raw_exons.values() if isinstance(v, dict)]

        if not tx_items:
            return None, None

        parsed_txs: list[tuple[int, int, dict[str, Any], int]] = []
        for idx, tx in enumerate(tx_items, start=1):
            raw_chr = str(tx.get("chr") or "").replace("chr", "")
            # Match chromosome key case-insensitively against species_cfg.nc_chrom_map
            matched_chr = next(
                (k for k in species_cfg.nc_chrom_map if k.upper() == raw_chr.upper()),
                None,
            )
            if not matched_chr:
                continue
            cds_s = tx.get("cdsstart")
            cds_e = tx.get("cdsend")
            pos_list = tx.get("position") or []
            if cds_s is None or cds_e is None or not pos_list:
                continue
            cds_s_1 = int(cds_s) + 1  # UCSC 0-based -> 1-based inclusive
            cds_e_1 = int(cds_e)
            if cds_s_1 >= cds_e_1:
                continue

            strand = int(tx.get("strand") or 1)
            acc = str(tx.get("transcript") or f"{gene_symbol}_tx{idx}")
            exons_out: list[dict[str, Any]] = []
            coding_bp = 0
            for e_i, pair in enumerate(pos_list, start=1):
                if not isinstance(pair, (list, tuple)) or len(pair) < 2:
                    continue
                es = int(pair[0]) + 1
                ee = int(pair[1])
                if es > ee:
                    es, ee = ee, es
                exons_out.append({"id": f"{acc}.E{e_i}", "start": es, "end": ee, "strand": strand})
                ovs, ove = max(es, cds_s_1), min(ee, cds_e_1)
                if ovs <= ove:
                    coding_bp += ove - ovs + 1

            if coding_bp < 9:
                continue
            prot_aa = max(1, (coding_bp // 3) - 1)  # RefSeq cds includes terminal stop codon
            tx_s = min(e["start"] for e in exons_out)
            tx_e = max(e["end"] for e in exons_out)
            txstart_0 = tx.get("txstart")
            txend_1 = tx.get("txend")
            if txstart_0 is not None and txend_1 is not None:
                tss_1 = (int(txstart_0) + 1) if strand >= 0 else int(txend_1)
            else:
                tss_1 = tx_s if strand >= 0 else tx_e

            is_nm = 0 if acc.upper().startswith("NM_") else 1
            tx_dict = {
                "id": acc,
                "object_type": "Transcript",
                "display_name": acc,
                "biotype": "protein_coding",
                "is_canonical": 0,
                "is_mane_select": False,
                "gencode_primary": 1 if acc.upper().startswith("NM_") else 0,
                "seq_region_name": matched_chr,
                "strand": strand,
                "start": tx_s,
                "end": tx_e,
                "Translation": {
                    "id": f"{acc}_P",
                    "start": cds_s_1,
                    "end": cds_e_1,
                    "length": prot_aa,
                },
                "Exon": exons_out,
            }
            parsed_txs.append((is_nm, -prot_aa, tx_dict, tss_1))

        if not parsed_txs:
            return None, None

        parsed_txs.sort(key=lambda x: (x[0], x[1], str(x[2]["id"])))
        primary_chr = str(parsed_txs[0][2]["seq_region_name"])
        primary_strand = int(parsed_txs[0][2]["strand"])
        primary_tss = parsed_txs[0][3]

        kept_txs: list[dict[str, Any]] = []
        for idx, (_nm, _neg_len, t_d, _tss) in enumerate(parsed_txs):
            if t_d["seq_region_name"] != primary_chr or t_d["strand"] != primary_strand:
                continue
            if idx == 0:
                t_d["is_canonical"] = 1
                t_d["is_mane_select"] = True
            kept_txs.append(t_d)

        gene_obj = {
            "id": gene_id,
            "object_type": "Gene",
            "assembly_name": species_cfg.assembly,
            "display_name": gene_symbol,
            "description": description,
            "seq_region_name": primary_chr,
            "strand": primary_strand,
            "start": min(int(t["start"]) for t in kept_txs),
            "end": max(int(t["end"]) for t in kept_txs),
            "Transcript": kept_txs,
        }
        return gene_obj, primary_tss

    def build_gene_from_ncbi_product_report(
        self,
        entrez_id: int,
        gene_id: str,
        gene_symbol: str,
        description: str,
        species_cfg: SpeciesConfig,
    ) -> tuple[Optional[dict[str, Any]], Optional[int]]:
        """Build an Ensembl-compatible Gene dict + TSS from NCBI Datasets v2 `product_report` (~300 ms)."""
        url = f"{NCBI_DATASETS_GENE_URL}/{int(entrez_id)}/product_report"
        try:
            data = self.http.fetch_json(url, max_retries=1, timeout_override=6.5)
        except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
            return None, None

        reports = (data or {}).get("reports") if isinstance(data, dict) else []
        if not reports or not isinstance(reports[0], dict):
            return None, None
        prod = reports[0].get("product") or {}
        raw_txs = prod.get("transcripts") or []
        nc_to_chr = {acc: ch for ch, acc in species_cfg.nc_chrom_map.items()}

        parsed_txs: list[tuple[int, int, int, dict[str, Any], int]] = []
        for tx in raw_txs:
            if not isinstance(tx, dict) or tx.get("type") != "PROTEIN_CODING":
                continue
            acc = str(tx.get("accession_version") or "").strip()
            cds_ranges = ((tx.get("cds") or {}).get("range")) or []
            gls = tx.get("genomic_locations") or []
            if not acc or not cds_ranges or not gls:
                continue
            cds_cdna_beg = int(cds_ranges[0]["begin"])
            cds_cdna_end = int(cds_ranges[0]["end"])

            matched_gl = None
            matched_chr = None
            for gl in gls:
                if not isinstance(gl, dict):
                    continue
                g_acc = str(gl.get("genomic_accession_version") or "")
                if g_acc in nc_to_chr and gl.get("exons"):
                    matched_gl = gl
                    matched_chr = nc_to_chr[g_acc]
                    break
            if not matched_gl or not matched_chr:
                continue

            gr = matched_gl.get("genomic_range") or {}
            strand = -1 if gr.get("orientation") == "minus" else 1
            raw_exons = sorted(
                matched_gl.get("exons") or [],
                key=lambda e: int(e.get("order") or 0),
            )
            if not raw_exons:
                continue

            # Map 1-based cDNA [cds_cdna_beg, cds_cdna_end] onto 1-based genomic coordinates across 5'->3' ordered exons
            cum = 0
            g_cds_points: list[int] = []
            exons_out: list[dict[str, Any]] = []
            for e_idx, ex in enumerate(raw_exons, start=1):
                eb = min(int(ex["begin"]), int(ex["end"]))
                ee = max(int(ex["begin"]), int(ex["end"]))
                elen = ee - eb + 1
                exons_out.append({"id": f"{acc}.E{e_idx}", "start": eb, "end": ee, "strand": strand})
                ex_c_start = cum + 1
                ex_c_end = cum + elen
                for target_cdna in (cds_cdna_beg, cds_cdna_end):
                    if ex_c_start <= target_cdna <= ex_c_end:
                        offset = target_cdna - ex_c_start
                        g_pos = (eb + offset) if strand >= 0 else (ee - offset)
                        g_cds_points.append(g_pos)
                cum += elen

            if len(g_cds_points) < 2:
                continue
            t_gstart, t_gend = min(g_cds_points), max(g_cds_points)
            prot_len = int(((tx.get("protein") or {}).get("length")) or max(1, (cds_cdna_end - cds_cdna_beg + 1) // 3 - 1))
            tx_s = min(e["start"] for e in exons_out)
            tx_e = max(e["end"] for e in exons_out)
            tss_1 = tx_s if strand >= 0 else tx_e

            sel_cat = str(tx.get("select_category") or "")
            is_sel = 0 if ("SELECT" in sel_cat or "MANE" in sel_cat) else 1
            is_nm = 0 if acc.upper().startswith("NM_") else 1

            tx_dict = {
                "id": acc.split(".")[0],
                "object_type": "Transcript",
                "display_name": acc,
                "biotype": "protein_coding",
                "is_canonical": 0,
                "is_mane_select": False,
                "gencode_primary": 1 if is_nm == 0 else 0,
                "seq_region_name": matched_chr,
                "strand": strand,
                "start": tx_s,
                "end": tx_e,
                "Translation": {
                    "id": str(((tx.get("protein") or {}).get("accession_version")) or f"{acc}_P"),
                    "start": t_gstart,
                    "end": t_gend,
                    "length": prot_len,
                },
                "Exon": exons_out,
            }
            parsed_txs.append((is_sel, is_nm, -prot_len, tx_dict, tss_1))

        if not parsed_txs:
            return None, None

        parsed_txs.sort(key=lambda x: (x[0], x[1], x[2], str(x[3]["display_name"])))
        primary_chr = str(parsed_txs[0][3]["seq_region_name"])
        primary_strand = int(parsed_txs[0][3]["strand"])
        primary_tss = parsed_txs[0][4]

        kept_txs: list[dict[str, Any]] = []
        for idx, (_sel, _nm, _nlen, t_d, _tss) in enumerate(parsed_txs):
            if t_d["seq_region_name"] != primary_chr or t_d["strand"] != primary_strand:
                continue
            if idx == 0:
                t_d["is_canonical"] = 1
                t_d["is_mane_select"] = True
            kept_txs.append(t_d)

        gene_obj = {
            "id": gene_id,
            "object_type": "Gene",
            "assembly_name": species_cfg.assembly,
            "display_name": gene_symbol,
            "description": description,
            "seq_region_name": primary_chr,
            "strand": primary_strand,
            "start": min(int(t["start"]) for t in kept_txs),
            "end": max(int(t["end"]) for t in kept_txs),
            "Transcript": kept_txs,
        }
        return gene_obj, primary_tss

    def fetch_alliance_paralogs_and_ortholog(
        self,
        mod_curie: Optional[str],
        min_paralog_identity: float = 35.0,
    ) -> tuple[list[dict[str, Any]], Optional[str]]:
        """Fetch within-species paralogs (with % identity) and stringent 1:1 Human Ortholog via Alliance (~150 ms)."""
        if not mod_curie:
            return [], None

        par_url = f"{ALLIANCE_GENE_URL}/{urllib.parse.quote(mod_curie)}/paralogs?limit=15"
        ort_url = f"{ALLIANCE_GENE_URL}/{urllib.parse.quote(mod_curie)}/orthologs?limit=35"

        def _get(u: str) -> dict[str, Any]:
            try:
                res = self.http.fetch_json(u, max_retries=0, timeout_override=4.0)
                return res if isinstance(res, dict) else {}
            except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
                return {}

        with ThreadPoolExecutor(max_workers=2) as pool:
            f_p = pool.submit(_get, par_url)
            f_o = pool.submit(_get, ort_url)
            p_data = f_p.result()
            o_data = f_o.result()

        paralogs: list[dict[str, Any]] = []
        seen_syms: set[str] = set()
        for row in p_data.get("results") or []:
            if not isinstance(row, dict):
                continue
            gtgp = row.get("geneToGeneParalogy") or {}
            obj_g = gtgp.get("objectGene") or {}
            sym = str(((obj_g.get("geneSymbol") or {}).get("displayText")) or "").strip()
            ext_id = str(obj_g.get("primaryExternalId") or sym).strip()
            ident = float(gtgp.get("identity") or gtgp.get("similarity") or row.get("identity") or 0.0)
            sim = float(gtgp.get("similarity") or row.get("similarity") or ident)
            if not sym or sym in seen_syms or ident < min_paralog_identity:
                continue
            seen_syms.add(sym)
            paralogs.append(
                {
                    "gene_symbol": sym,
                    "ensembl_id": ext_id,
                    "percent_identity": round(ident, 1),
                    "target_percent_identity": round(sim, 1),
                    "ancestor_taxonomy": "Alliance DIOPT",
                }
            )

        # Find best stringent human ortholog (NCBITaxon:9606)
        human_candidates: list[tuple[int, int, str]] = []
        for row in o_data.get("results") or []:
            if not isinstance(row, dict):
                continue
            orth = row.get("geneToGeneOrthologyGenerated") or {}
            obj_g = orth.get("objectGene") or {}
            taxon_curie = str(((obj_g.get("taxon") or {}).get("curie")) or "")
            if taxon_curie != "NCBITaxon:9606":
                continue
            h_sym = str(((obj_g.get("geneSymbol") or {}).get("displayText")) or "").strip()
            if not h_sym:
                continue
            is_best = 0 if ((orth.get("isBestScore") or {}).get("name") == "Yes") else 1
            n_methods = len(orth.get("predictionMethodsMatched") or [])
            human_candidates.append((is_best, -n_methods, h_sym))

        human_orth_sym: Optional[str] = None
        if human_candidates:
            human_candidates.sort()
            human_orth_sym = human_candidates[0][2]

        return paralogs, human_orth_sym
