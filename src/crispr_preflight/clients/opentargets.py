"""Open Targets Platform GraphQL v4 client for DepMap essentiality, chemical probes, and drugs."""

from __future__ import annotations

import statistics
import urllib.error
from typing import Any, Optional

from crispr_preflight.cache import CachedHttpClient, NetworkFetchError
from crispr_preflight.models import ChemicalTool, EssentialityProfile

OPENTARGETS_URL = "https://api.platform.opentargets.org/api/v4/graphql"


class OpenTargetsClient:
    """Queries Open Targets v4 for cell-line Chronos scores and orthogonal chemical tools."""

    def __init__(self, http: CachedHttpClient) -> None:
        self.http = http

    def _graphql(self, query: str, variables: dict[str, Any]) -> Optional[dict[str, Any]]:
        try:
            resp = self.http.fetch_json(
                OPENTARGETS_URL,
                method="POST",
                json_body={"query": query, "variables": variables},
                max_retries=1,
                timeout_override=9.0,
            )
        except (urllib.error.HTTPError, urllib.error.URLError, NetworkFetchError):
            return None
        if not isinstance(resp, dict) or "errors" in resp:
            return None
        return resp.get("data")

    @staticmethod
    def _extract_swissprot_id(protein_ids: Any) -> Optional[str]:
        for item in protein_ids or []:
            if isinstance(item, dict) and item.get("source") == "uniprot_swissprot" and item.get("id"):
                return str(item["id"]).split("-")[0]
        return None

    def resolve_target_metadata(self, gene_query: str) -> Optional[dict[str, Optional[str]]]:
        """Fast (~100-250 ms) symbol/alias/ENSG -> metadata resolution via Open Targets."""
        cleaned = gene_query.strip()
        if not cleaned:
            return None
        if cleaned.upper().startswith("ENSG"):
            base_id = cleaned.upper().split(".")[0]
            q_id = """
            query ResolveById($ensemblId: String!) {
              target(ensemblId: $ensemblId) {
                id
                approvedSymbol
                approvedName
                proteinIds {
                  id
                  source
                }
              }
            }
            """
            data = self._graphql(q_id, {"ensemblId": base_id})
            tgt = (data or {}).get("target")
            if isinstance(tgt, dict) and tgt.get("id"):
                return {
                    "ensembl_id": str(tgt["id"]),
                    "gene_symbol": str(tgt.get("approvedSymbol") or base_id),
                    "approved_name": str(tgt.get("approvedName") or ""),
                    "swissprot_id": self._extract_swissprot_id(tgt.get("proteinIds")),
                }
            return {
                "ensembl_id": base_id,
                "gene_symbol": base_id,
                "approved_name": "",
                "swissprot_id": None,
            }

        q = """
        query ResolveGene($q: String!) {
          search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 10}) {
            hits {
              id
              name
              object {
                ... on Target {
                  id
                  approvedSymbol
                  approvedName
                  proteinIds {
                    id
                    source
                  }
                  symbolSynonyms {
                    label
                  }
                }
              }
            }
          }
        }
        """
        data = self._graphql(q, {"q": cleaned})
        hits = (((data or {}).get("search") or {}).get("hits")) or []
        if not hits:
            return None

        target_upper = cleaned.upper()
        # 1. Exact HGNC approvedSymbol match
        for h in hits:
            obj = h.get("object") or {}
            ensg = str(obj.get("id") or h.get("id") or "")
            sym = str(obj.get("approvedSymbol") or h.get("name") or "")
            if ensg.startswith("ENSG") and sym.upper() == target_upper:
                return {
                    "ensembl_id": ensg,
                    "gene_symbol": sym,
                    "approved_name": str(obj.get("approvedName") or ""),
                    "swissprot_id": self._extract_swissprot_id(obj.get("proteinIds")),
                }

        # 2. Exact synonym match (e.g. BRG1 -> SMARCA4, HER2 -> ERBB2, P53 -> TP53)
        for h in hits:
            obj = h.get("object") or {}
            ensg = str(obj.get("id") or h.get("id") or "")
            sym = str(obj.get("approvedSymbol") or h.get("name") or "")
            syns = [
                str(s.get("label") or "").upper()
                for s in (obj.get("symbolSynonyms") or [])
                if isinstance(s, dict)
            ]
            if ensg.startswith("ENSG") and target_upper in syns:
                return {
                    "ensembl_id": ensg,
                    "gene_symbol": sym,
                    "approved_name": str(obj.get("approvedName") or ""),
                    "swissprot_id": self._extract_swissprot_id(obj.get("proteinIds")),
                }

        return None

    def fetch_target_bundle(
        self,
        ensembl_id: str,
        cell_line_query: Optional[str] = None,
        tissue_query: Optional[str] = None,
        max_tools: int = 6,
        min_paralog_identity: float = 35.0,
    ) -> tuple[EssentialityProfile, Optional[dict[str, float]], list[ChemicalTool], Optional[list[dict[str, Any]]]]:
        """Single-shot GraphQL fetch combining DepMap Chronos, gnomAD constraint, paralogs, probes, and drugs."""
        bundle_query = """
        query TargetPreflightBundle($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            id
            isEssential
            geneticConstraint {
              constraintType
              oe
              oeUpper
              score
            }
            homologues {
              speciesId
              homologyType
              targetGeneId
              targetGeneSymbol
              queryPercentageIdentity
              targetPercentageIdentity
            }
            depMapEssentiality {
              tissueName
              screens {
                depmapId
                cellLineName
                diseaseFromSource
                geneEffect
              }
            }
            chemicalProbes {
              id
              control
              drugId
              mechanismOfAction
              isHighQuality
              probesDrugsScore
              scoreInCells
              urls {
                niceName
                url
              }
            }
            drugAndClinicalCandidates {
              count
              rows {
                maxClinicalStage
                drug {
                  id
                  name
                  drugType
                  mechanismsOfAction {
                    rows {
                      mechanismOfAction
                    }
                  }
                }
              }
            }
          }
        }
        """
        data = self._graphql(bundle_query, {"ensemblId": ensembl_id})
        if not data or not data.get("target"):
            ess, const = self.fetch_essentiality(ensembl_id, cell_line_query, tissue_query)
            tools = self.fetch_chemical_tools(ensembl_id, max_tools=max_tools)
            return ess, const, tools, None

        target = data["target"]
        ess, const = self._parse_essentiality_target(target, cell_line_query, tissue_query)
        tools = self._parse_chemical_tools_target(target, max_tools=max_tools)
        paralogs = self._parse_human_paralogs_target(target, ensembl_id, min_identity=min_paralog_identity)
        return ess, const, tools, paralogs

    @staticmethod
    def _parse_human_paralogs_target(
        target: dict[str, Any],
        ensembl_id: str,
        min_identity: float = 35.0,
    ) -> Optional[list[dict[str, Any]]]:
        if "homologues" not in target or target["homologues"] is None:
            return None
        raw_out: list[dict[str, Any]] = []
        for h in target.get("homologues") or []:
            if str(h.get("speciesId") or "") != "9606":
                continue
            htype = str(h.get("homologyType") or "").lower()
            if "paralog" not in htype:
                continue
            tid = str(h.get("targetGeneId") or "")
            if not tid or tid == ensembl_id:
                continue
            q_id = float(h.get("queryPercentageIdentity") or 0.0)
            t_id = float(h.get("targetPercentageIdentity") or 0.0)
            if max(q_id, t_id) < min_identity:
                continue
            sym = str(h.get("targetGeneSymbol") or tid)
            raw_out.append(
                {
                    "ensembl_id": tid,
                    "gene_symbol": sym,
                    "percent_identity": round(q_id, 1),
                    "target_percent_identity": round(t_id, 1),
                    "ancestor_taxonomy": h.get("homologyType"),
                }
            )
        dedup: dict[str, dict[str, Any]] = {}
        for item in sorted(raw_out, key=lambda x: x["percent_identity"], reverse=True):
            if item["gene_symbol"] not in dedup:
                dedup[item["gene_symbol"]] = item
        return list(dedup.values())

    def _parse_essentiality_target(
        self,
        target: dict[str, Any],
        cell_line_query: Optional[str] = None,
        tissue_query: Optional[str] = None,
    ) -> tuple[EssentialityProfile, Optional[dict[str, float]]]:
        profile = EssentialityProfile()
        ot_constraint: dict[str, float] = {}

        if target.get("isEssential") is not None:
            profile.is_core_essential_ot = bool(target.get("isEssential"))

        for gc in target.get("geneticConstraint") or []:
            ctype = str(gc.get("constraintType", "")).lower()
            if ctype == "lof":
                if gc.get("oe") is not None:
                    ot_constraint["oe_lof"] = float(gc["oe"])
                if gc.get("oeUpper") is not None:
                    ot_constraint["loeuf"] = float(gc["oeUpper"])
                if gc.get("score") is not None:
                    ot_constraint["pli"] = float(gc["score"])

        all_effects: list[float] = []
        lineage_effects: list[float] = []
        matched_lineage_label: Optional[str] = None
        best_cell_match: Optional[tuple[int, str, str, float]] = None

        norm_cell = self._norm_cell_name(cell_line_query) if cell_line_query else ""
        norm_tissue = tissue_query.strip().lower() if tissue_query else ""

        for tissue_group in target.get("depMapEssentiality") or []:
            t_name = str(tissue_group.get("tissueName") or "Unknown")
            is_tissue_hit = bool(
                norm_tissue
                and (
                    norm_tissue in t_name.lower()
                    or t_name.lower() in norm_tissue
                )
            )
            if is_tissue_hit and matched_lineage_label is None:
                matched_lineage_label = t_name

            for scr in tissue_group.get("screens") or []:
                eff = scr.get("geneEffect")
                if eff is None:
                    continue
                eff_val = float(eff)
                all_effects.append(eff_val)
                if is_tissue_hit:
                    lineage_effects.append(eff_val)

                if norm_cell:
                    cname = str(scr.get("cellLineName") or "")
                    did = str(scr.get("depmapId") or "")
                    norm_cname = self._norm_cell_name(cname)
                    if norm_cname == norm_cell or did.lower() == norm_cell:
                        cand = (0, cname or did, t_name, eff_val)
                        if best_cell_match is None or cand[0] < best_cell_match[0]:
                            best_cell_match = cand
                    elif norm_cell in norm_cname or (len(norm_cname) >= 3 and norm_cname in norm_cell):
                        cand = (1, cname or did, t_name, eff_val)
                        if best_cell_match is None or cand[0] < best_cell_match[0]:
                            best_cell_match = cand

        if all_effects:
            profile.depmap_total_lines = len(all_effects)
            profile.depmap_essential_lines = sum(1 for v in all_effects if v <= -1.0)
            profile.depmap_pan_essential_fraction = round(
                profile.depmap_essential_lines / len(all_effects), 3
            )
            profile.depmap_median_chronos = round(statistics.median(all_effects), 3)

        if best_cell_match is not None:
            _, cname, t_name, eff_val = best_cell_match
            profile.matched_cell_line = cname
            profile.matched_cell_line_lineage = t_name
            profile.matched_cell_line_chronos = round(eff_val, 3)

        if lineage_effects and matched_lineage_label:
            profile.matched_lineage = matched_lineage_label
            profile.matched_lineage_median_chronos = round(statistics.median(lineage_effects), 3)
            profile.matched_lineage_cell_count = len(lineage_effects)

        return profile, (ot_constraint or None)

    def _parse_chemical_tools_target(
        self,
        target: dict[str, Any],
        max_tools: int = 8,
    ) -> list[ChemicalTool]:
        tools: list[ChemicalTool] = []
        seen_names: set[str] = set()

        probes = target.get("chemicalProbes") or []
        probes_sorted = sorted(
            probes,
            key=lambda p: (
                0 if p.get("isHighQuality") else 1,
                -(float(p.get("probesDrugsScore") or p.get("scoreInCells") or 0.0)),
            ),
        )
        for p in probes_sorted:
            pname = str(p.get("id") or p.get("drugId") or "").strip()
            if not pname or pname.lower() in seen_names:
                continue
            seen_names.add(pname.lower())
            moa_raw = p.get("mechanismOfAction")
            if isinstance(moa_raw, list):
                moa = ", ".join(str(x) for x in moa_raw if x) or "Chemical probe modulator"
            else:
                moa = str(moa_raw or "Chemical probe modulator")
            hq = bool(p.get("isHighQuality"))
            cell_score = p.get("scoreInCells")
            qual_parts = ["High-quality probe" if hq else "Candidate probe"]
            if cell_score is not None:
                qual_parts.append(f"cell score: {cell_score}")
            has_ctrl = bool(p.get("control"))
            url = None
            for u in p.get("urls") or []:
                if isinstance(u, dict) and u.get("url"):
                    url = str(u["url"])
                    break
            tools.append(
                ChemicalTool(
                    name=pname,
                    identifier=str(p.get("drugId") or pname),
                    tool_type="Chemical Probe",
                    mechanism_of_action=moa,
                    quality_or_stage=" | ".join(qual_parts),
                    has_negative_control=has_ctrl,
                    source_url=url,
                )
            )
            if len(tools) >= max_tools:
                return tools

        stage_rank = {"APPROVAL": 5, "PHASE_4": 5, "PHASE_3": 4, "PHASE_2": 3, "PHASE_1": 2, "PHASE_0": 1}
        cand_obj = target.get("drugAndClinicalCandidates") or {}
        rows = cand_obj.get("rows") or []
        rows_sorted = sorted(
            rows,
            key=lambda r: -stage_rank.get(str(r.get("maxClinicalStage") or ""), 0),
        )
        for row in rows_sorted:
            drug = row.get("drug") or {}
            dname = str(drug.get("name") or "").strip()
            did = str(drug.get("id") or "").strip()
            if not dname or dname.lower() in seen_names:
                continue
            seen_names.add(dname.lower())
            raw_stage = str(row.get("maxClinicalStage") or "Clinical")
            stage = "Approved" if raw_stage == "APPROVAL" else raw_stage.replace("_", " ").title()
            dtype = str(drug.get("drugType") or "Small molecule")
            moa_rows = ((drug.get("mechanismsOfAction") or {}).get("rows")) or []
            moa = str(moa_rows[0].get("mechanismOfAction")) if moa_rows else dtype
            tools.append(
                ChemicalTool(
                    name=dname,
                    identifier=did,
                    tool_type=f"Clinical ({dtype})",
                    mechanism_of_action=moa,
                    quality_or_stage=stage,
                    has_negative_control=False,
                    source_url=f"https://platform.opentargets.org/drug/{did}" if did else None,
                )
            )
            if len(tools) >= max_tools:
                break

        return tools

    def fetch_essentiality(
        self,
        ensembl_id: str,
        cell_line_query: Optional[str] = None,
        tissue_query: Optional[str] = None,
    ) -> tuple[EssentialityProfile, Optional[dict[str, float]]]:
        """Fetch DepMap CRISPR Chronos scores and fallback constraint metrics from Open Targets."""
        query = """
        query TargetEssentiality($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            id
            isEssential
            geneticConstraint {
              constraintType
              oe
              oeUpper
              score
            }
            depMapEssentiality {
              tissueName
              screens {
                depmapId
                cellLineName
                diseaseFromSource
                geneEffect
              }
            }
          }
        }
        """
        data = self._graphql(query, {"ensemblId": ensembl_id})
        if not data or not data.get("target"):
            return EssentialityProfile(), None
        return self._parse_essentiality_target(data["target"], cell_line_query, tissue_query)

    def fetch_chemical_tools(self, ensembl_id: str, max_tools: int = 8) -> list[ChemicalTool]:
        """Fetch curated chemical probes and clinical small-molecule/degrader modulators."""
        query = """
        query TargetChemicalTools($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            chemicalProbes {
              id
              control
              drugId
              mechanismOfAction
              isHighQuality
              probesDrugsScore
              scoreInCells
              urls {
                niceName
                url
              }
            }
            drugAndClinicalCandidates {
              count
              rows {
                maxClinicalStage
                drug {
                  id
                  name
                  drugType
                  mechanismsOfAction {
                    rows {
                      mechanismOfAction
                    }
                  }
                }
              }
            }
          }
        }
        """
        data = self._graphql(query, {"ensemblId": ensembl_id})
        if not data or not data.get("target"):
            return []
        return self._parse_chemical_tools_target(data["target"], max_tools=max_tools)

    @staticmethod
    def _norm_cell_name(name: str) -> str:
        """Strip hyphens, spaces, and punctuation so 'HEK-293T' matches 'HEK293T'."""
        return "".join(ch for ch in name.lower() if ch.isalnum())
