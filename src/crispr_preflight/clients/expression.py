"""GTEx Portal v2 and Human Protein Atlas (HPA) clients for tissue expression & localization."""

from __future__ import annotations

import urllib.error
import urllib.parse
from typing import Any, Optional

from crispr_preflight.cache import CachedHttpClient, NetworkFetchError

GTEX_BASE = "https://gtexportal.org/api/v2"
GTEX_DATASET = "gtex_v10"
GTEX_GENCODE = "v39"
HPA_BASE = "https://www.proteinatlas.org"

# Convenient bench shorthand -> GTEx tissueSiteDetailId substring hints
TISSUE_ALIASES: dict[str, list[str]] = {
    "blood": ["whole_blood", "cells_ebv-transformed_lymphocytes"],
    "pbmc": ["whole_blood"],
    "myeloid": ["whole_blood", "spleen"],
    "lymphoid": ["cells_ebv-transformed_lymphocytes", "spleen"],
    "lung": ["lung"],
    "liver": ["liver"],
    "kidney": ["kidney_cortex", "kidney_medulla"],
    "brain": ["brain_cortex", "brain_frontal_cortex_ba9", "brain_cerebellum"],
    "colon": ["colon_sigmoid", "colon_transverse"],
    "gut": ["colon_sigmoid", "small_intestine_terminal_ileum"],
    "intestine": ["small_intestine_terminal_ileum", "colon_transverse"],
    "breast": ["breast_mammary_tissue"],
    "skin": ["skin_not_sun_exposed_suprapubic", "skin_sun_exposed_lower_leg"],
    "pancreas": ["pancreas"],
    "prostate": ["prostate"],
    "ovary": ["ovary"],
    "stomach": ["stomach"],
    "heart": ["heart_left_ventricle", "heart_atrial_appendage"],
    "muscle": ["muscle_skeletal"],
    "adipose": ["adipose_subcutaneous", "adipose_visceral_omentum"],
    "fibroblast": ["cells_cultured_fibroblasts"],
}


class ExpressionClient:
    """Queries GTEx v2 for quantitative mRNA TPM and HPA for subcellular localization & IHC."""

    def __init__(self, http: CachedHttpClient) -> None:
        self.http = http

    def fetch_gtex_gene_medians(self, gene_symbol: str) -> dict[str, float]:
        """Return {tissueSiteDetailId: median_tpm} across GTEx v10 tissues for a gene symbol."""
        ref_url = (
            f"{GTEX_BASE}/reference/gene?"
            f"{urllib.parse.urlencode({'geneId': gene_symbol, 'gencodeVersion': GTEX_GENCODE})}"
        )
        try:
            ref_resp = self.http.fetch_json(ref_url, max_retries=1, timeout_override=7.0) or {}
        except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
            return {}

        genes = ref_resp.get("data") if isinstance(ref_resp, dict) else ref_resp
        if not genes or not isinstance(genes, list):
            return {}

        best = genes[0]
        for g in genes:
            if str(g.get("geneSymbol", "")).lower() == gene_symbol.lower():
                best = g
                break

        gencode_id = best.get("gencodeId")
        if not gencode_id:
            return {}

        expr_url = (
            f"{GTEX_BASE}/expression/medianGeneExpression?"
            f"{urllib.parse.urlencode({'gencodeId': gencode_id, 'datasetId': GTEX_DATASET, 'itemsPerPage': 250})}"
        )
        try:
            expr_resp = self.http.fetch_json(expr_url, max_retries=1, timeout_override=7.0) or {}
        except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
            return {}

        rows = expr_resp.get("data") if isinstance(expr_resp, dict) else expr_resp
        if not rows or not isinstance(rows, list):
            return {}

        out: dict[str, float] = {}
        for row in rows:
            tid = str(row.get("tissueSiteDetailId") or "")
            med = row.get("median")
            if tid and med is not None:
                out[tid] = round(float(med), 2)
        return out

    @staticmethod
    def match_tissue_tpm(
        medians: dict[str, float],
        tissue_query: Optional[str],
    ) -> tuple[Optional[str], Optional[float]]:
        """Match a user tissue string (e.g. 'lung', 'blood', 'Brain_Cortex') to GTEx tissue ID."""
        if not tissue_query or not medians:
            return None, None
        q = tissue_query.strip().lower()
        q_norm = q.replace(" ", "_").replace("-", "_")

        # 1. Exact ID match
        for tid, val in medians.items():
            if tid.lower() == q_norm:
                return tid, val

        # 2. Alias lookup
        for alias_key, candidates in TISSUE_ALIASES.items():
            if q == alias_key or alias_key in q:
                for cand in candidates:
                    for tid, val in medians.items():
                        if tid.lower() == cand:
                            return tid, val

        # 3. Substring match (pick highest-expressed matching sub-tissue if several match, e.g. brain)
        sub_hits = [
            (tid, val)
            for tid, val in medians.items()
            if q_norm in tid.lower() or q in tid.lower().replace("_", " ")
        ]
        if sub_hits:
            sub_hits.sort(key=lambda x: x[1], reverse=True)
            return sub_hits[0]

        return None, None

    def fetch_hpa_annotations(self, ensembl_id: str) -> dict[str, Any]:
        """Fetch subcellular localization and RNA tissue specificity from the Human Protein Atlas."""
        result: dict[str, Any] = {
            "subcellular_locations": [],
            "rna_specificity": None,
            "ihc_level": None,
        }
        params = urllib.parse.urlencode(
            {
                "search": ensembl_id,
                "columns": "g,eg,scl,scml,rnats",
                "format": "json",
                "compress": "no",
            }
        )
        search_url = f"{HPA_BASE}/api/search_download.php?{params}"
        try:
            rows = self.http.fetch_json(search_url, max_retries=1, timeout_override=6.0) or []
            if isinstance(rows, list) and rows:
                matched = rows[0]
                for r in rows:
                    if str(r.get("Ensembl", "")) == ensembl_id:
                        matched = r
                        break
                locs = matched.get("Subcellular main location") or matched.get("Subcellular location") or []
                if isinstance(locs, str):
                    locs = [s.strip() for s in locs.split(",") if s.strip()]
                result["subcellular_locations"] = list(locs)
                result["rna_specificity"] = matched.get("RNA tissue specificity")
        except (urllib.error.HTTPError, urllib.error.URLError, ValueError, NetworkFetchError):
            pass

        return result
