"""gnomAD v4 GraphQL client for human loss-of-function constraint (pLI and LOEUF)."""

from __future__ import annotations

import urllib.error
from typing import Optional

from crispr_preflight.cache import CachedHttpClient, NetworkFetchError

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"


class GnomadClient:
    """Fetches gene-level loss-of-function constraint metrics from gnomAD GRCh38."""

    def __init__(self, http: CachedHttpClient) -> None:
        self.http = http

    def fetch_constraint(self, gene_symbol: str) -> Optional[dict[str, Optional[float]]]:
        """Return {'pli': ..., 'loeuf': ..., 'oe_lof': ...} or None if unavailable."""
        query = """
        query GeneConstraint($geneSymbol: String!, $referenceGenome: ReferenceGenomeId!) {
          gene(gene_symbol: $geneSymbol, reference_genome: $referenceGenome) {
            gene_id
            symbol
            gnomad_constraint {
              pli
              oe_lof
              oe_lof_upper
              oe_mis
            }
          }
        }
        """
        try:
            resp = self.http.fetch_json(
                GNOMAD_API_URL,
                method="POST",
                json_body={
                    "query": query,
                    "variables": {"geneSymbol": gene_symbol, "referenceGenome": "GRCh38"},
                },
                max_retries=0,
                timeout_override=6.5,
            )
        except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
            return None

        if not isinstance(resp, dict):
            return None
        gene_data = (resp.get("data") or {}).get("gene")
        if not gene_data or not gene_data.get("gnomad_constraint"):
            return None

        c = gene_data["gnomad_constraint"]
        return {
            "pli": float(c["pli"]) if c.get("pli") is not None else None,
            "oe_lof": float(c["oe_lof"]) if c.get("oe_lof") is not None else None,
            "loeuf": float(c["oe_lof_upper"]) if c.get("oe_lof_upper") is not None else None,
        }
