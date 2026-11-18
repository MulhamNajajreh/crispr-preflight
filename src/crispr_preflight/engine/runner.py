"""Main pipeline orchestrator assembling a complete PreflightReport for a gene."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from typing import Any, Optional

from crispr_preflight.cache import CachedHttpClient
from crispr_preflight.clients.cellline_variants import CellLineAndVariantClient
from crispr_preflight.clients.ensembl import (
    GRCH38_NC_MAP,
    EnsemblClient,
)
from crispr_preflight.clients.expression import ExpressionClient
from crispr_preflight.clients.gnomad import GnomadClient
from crispr_preflight.clients.opentargets import OpenTargetsClient
from crispr_preflight.engine.exon_scorer import (
    design_crispri_tss_guides,
    reverse_complement,
    score_transcript_exons,
    select_coding_isoforms,
    select_top_gene_guides,
)
from crispr_preflight.engine.modality import classify_paralog_risks, evaluate_modality
from crispr_preflight.engine.validation import build_validation_plan
from crispr_preflight.models import (
    CellLineAlterationProfile,
    ExpressionProfile,
    GuideCandidate,
    PreflightReport,
    ProteinDomain,
)


class PreflightRunner:
    """Coordinates database lookups and biological rule engines to build a PreflightReport."""

    def __init__(
        self,
        http: Optional[CachedHttpClient] = None,
        check_pams: bool = True,
        query_gnomad_direct: bool = True,
    ) -> None:
        self.http = http or CachedHttpClient()
        self.ensembl = EnsemblClient(self.http)
        self.opentargets = OpenTargetsClient(self.http)
        self.gnomad = GnomadClient(self.http)
        self.expression = ExpressionClient(self.http)
        self.cellline_variants = CellLineAndVariantClient(self.http)
        self.check_pams = check_pams
        self.query_gnomad_direct = query_gnomad_direct

    def run(
        self,
        gene_query: str,
        tissue: Optional[str] = None,
        cell_line: Optional[str] = None,
    ) -> PreflightReport:
        """Execute end-to-end biological flight check for a human gene."""
        cleaned_query = gene_query.strip()
        if not cleaned_query:
            raise ValueError("Please provide a human gene symbol, alias, or Ensembl ID.")
        if len(cleaned_query) > 80:
            raise ValueError("Gene query is too long; provide a single human gene symbol or Ensembl ID.")

        sources_used: list[str] = []
        used_ensembl_rest = False

        # 1. Fast (~120-250 ms) symbol/alias -> (ENSG, symbol, UniProt Swiss-Prot ID) resolution via Open Targets
        meta = self.opentargets.resolve_target_metadata(cleaned_query)
        ebi_coord: Optional[dict[str, Any]] = None
        used_ebi = False
        preloaded_pfam_domains: list[ProteinDomain] = []
        mygene_tss: Optional[int] = None
        ccle_profile: Optional[CellLineAlterationProfile] = None

        if meta is not None and meta.get("ensembl_id"):
            fast_ensg = str(meta["ensembl_id"])
            fast_sym = str(meta.get("gene_symbol") or fast_ensg)
            approved_name = str(meta.get("approved_name") or "")
            swissprot_id = str(meta.get("swissprot_id") or "")

            def _resolve_fast_gene() -> tuple[dict[str, Any], Optional[dict[str, Any]], bool]:
                coord = self.ensembl.fetch_ebi_coordinates(swissprot_id) if swissprot_id else None
                g_data, u_ebi = self.ensembl.lookup_gene_fast(
                    fast_ensg,
                    fast_sym,
                    approved_name=approved_name,
                    ebi_coord=coord,
                )
                return g_data, coord, u_ebi

            with ThreadPoolExecutor(max_workers=8) as wave1:
                fut_gene = wave1.submit(_resolve_fast_gene)
                fut_ip = (
                    wave1.submit(self.ensembl.fetch_interpro_pfam_domains, swissprot_id, None)
                    if swissprot_id
                    else None
                )
                fut_ot = wave1.submit(
                    self.opentargets.fetch_target_bundle, fast_ensg, cell_line, tissue, 6, 35.0
                )
                fut_gtex = (
                    wave1.submit(self.expression.fetch_gtex_gene_medians, fast_sym)
                    if not fast_sym.startswith("ENSG")
                    else None
                )
                fut_hpa = wave1.submit(self.expression.fetch_hpa_annotations, fast_ensg)
                fut_tss = (
                    wave1.submit(self.cellline_variants.fetch_tss_from_mygene, fast_ensg)
                    if self.check_pams
                    else None
                )
                fut_ccle = (
                    wave1.submit(self.cellline_variants.fetch_ccle_alterations, fast_sym, cell_line)
                    if (cell_line and not fast_sym.startswith("ENSG"))
                    else None
                )

                gene_data, ebi_coord, used_ebi = fut_gene.result()
                if not used_ebi:
                    used_ensembl_rest = True
                if fut_ip is not None:
                    ip_doms = fut_ip.result()
                    preloaded_pfam_domains = self.ensembl.merge_uniprot_domain_features(
                        ip_doms, ebi_coord
                    )
                essentiality, ot_constraint, chem_tools, ot_paralogs = fut_ot.result()
                target_medians = fut_gtex.result() if fut_gtex is not None else None
                hpa_info = fut_hpa.result()
                mygene_tss = fut_tss.result() if fut_tss is not None else None
                ccle_profile = fut_ccle.result() if fut_ccle is not None else None
        else:
            gene_data = self.ensembl.lookup_gene_expanded(cleaned_query)
            used_ensembl_rest = True
            fast_ensg = str(gene_data["id"])
            fast_sym = str(gene_data.get("display_name") or cleaned_query.upper())
            swissprot_id = ""
            with ThreadPoolExecutor(max_workers=6) as wave1b:
                fut_ot = wave1b.submit(
                    self.opentargets.fetch_target_bundle, fast_ensg, cell_line, tissue, 6, 35.0
                )
                fut_gtex = wave1b.submit(self.expression.fetch_gtex_gene_medians, fast_sym)
                fut_hpa = wave1b.submit(self.expression.fetch_hpa_annotations, fast_ensg)
                fut_tss = (
                    wave1b.submit(self.cellline_variants.fetch_tss_from_mygene, fast_ensg)
                    if self.check_pams
                    else None
                )
                fut_ccle = (
                    wave1b.submit(self.cellline_variants.fetch_ccle_alterations, fast_sym, cell_line)
                    if cell_line
                    else None
                )
                essentiality, ot_constraint, chem_tools, ot_paralogs = fut_ot.result()
                target_medians = fut_gtex.result()
                hpa_info = fut_hpa.result()
                mygene_tss = fut_tss.result() if fut_tss is not None else None
                ccle_profile = fut_ccle.result() if fut_ccle is not None else None

        if used_ebi:
            sources_used.append("EMBL-EBI Proteins & InterPro APIs (GRCh38 / Pfam)")

        ensembl_id = str(gene_data["id"])
        gene_symbol = str(gene_data.get("display_name") or fast_sym or gene_query.upper())
        raw_desc = str(gene_data.get("description") or "").split("[Source:")[0].strip()
        chrom = str(gene_data.get("seq_region_name") or "")
        strand = int(gene_data.get("strand") or 1)

        if target_medians is None:
            target_medians = self.expression.fetch_gtex_gene_medians(gene_symbol)

        primary_raw, primary_summary, alt_raws, alt_summaries = select_coding_isoforms(gene_data)

        # Check if Swiss-Prot canonical sequence length matches primary isoform protein length (>98% of human genes)
        ebi_seq_len = len(str((ebi_coord or {}).get("sequence") or ""))
        can_use_fast_domains = bool(
            swissprot_id
            and ebi_seq_len > 0
            and ebi_seq_len == primary_summary.protein_length_aa
        )

        # 2. Prepare ordered 5'->3' coding intervals
        raw_exons = list(primary_raw.get("Exon") or [])
        ordered_exons = (
            sorted(raw_exons, key=lambda e: int(e["start"]))
            if strand >= 0
            else sorted(raw_exons, key=lambda e: int(e["end"]), reverse=True)
        )
        trans = primary_raw["Translation"]
        t_s, t_e = min(int(trans["start"]), int(trans["end"])), max(int(trans["start"]), int(trans["end"]))
        coding_intervals: list[tuple[int, int]] = []
        for ex in ordered_exons:
            es, ee = min(int(ex["start"]), int(ex["end"])), max(int(ex["start"]), int(ex["end"]))
            ovs, ove = max(es, t_s), min(ee, t_e)
            if ovs <= ove:
                coding_intervals.append((ovs, ove))

        # Resolve 5' Transcription Start Site (MyGene RefSeq TSS, or 5'-most primary exon coordinate)
        fallback_tss: Optional[int] = None
        if ordered_exons:
            fallback_tss = (
                min(int(ordered_exons[0]["start"]), int(ordered_exons[0]["end"]))
                if strand >= 0
                else max(int(ordered_exons[0]["start"]), int(ordered_exons[0]["end"]))
            )
        tss_pos = mygene_tss or fallback_tss

        flank_pad = 420  # +-420 bp around exons enables instant TIDE/ICE genotyping PCR primer design
        tss_iv: Optional[tuple[int, int]] = (
            (max(1, tss_pos - 60), tss_pos + 60) if (self.check_pams and tss_pos) else None
        )
        fetch_intervals = list(coding_intervals)
        if tss_iv is not None and tss_iv not in fetch_intervals:
            fetch_intervals.append(tss_iv)

        # Determine high-identity paralogs early from Open Targets bundle so paralog GTEx runs in Wave 2!
        resolved_tid, req_tpm = self.expression.match_tissue_tpm(target_medians, tissue)
        top_sorted = sorted(target_medians.items(), key=lambda x: x[1], reverse=True)[:5]
        ref_gtex_tissue = resolved_tid or (top_sorted[0][0] if top_sorted else None)

        high_id_paralog_syms: list[str] = []
        if ot_paralogs is not None:
            high_id_paralog_syms = [
                str(p["gene_symbol"])
                for p in ot_paralogs[:3]
                if max(float(p["percent_identity"]), float(p["target_percent_identity"])) >= 50.0
            ][:2]

        # 3. Wave 2: Run NCBI GRCh38 flanked exon + TSS sequences (+ fallback domains/CCLE/paralog GTEx) concurrently
        with ThreadPoolExecutor(max_workers=6) as wave2:
            fut_doms = (
                wave2.submit(self.ensembl.fetch_protein_domains, primary_summary.protein_id)
                if (not can_use_fast_domains and primary_summary.protein_id)
                else None
            )
            fut_flanks = (
                wave2.submit(self.ensembl.fetch_flanked_regions, chrom, fetch_intervals, flank_pad)
                if self.check_pams and fetch_intervals
                else None
            )
            fut_cds_fallback = (
                wave2.submit(self.ensembl.fetch_cds_sequence, primary_summary.transcript_id)
                if not self.check_pams
                else None
            )
            fut_ens_paralogs = (
                wave2.submit(self.ensembl.fetch_human_paralogs, ensembl_id, 35.0)
                if ot_paralogs is None
                else None
            )
            fut_ccle_late = (
                wave2.submit(self.cellline_variants.fetch_ccle_alterations, gene_symbol, cell_line)
                if (cell_line and ccle_profile is None)
                else None
            )
            paralog_gtex_futs = {
                sym: wave2.submit(self.expression.fetch_gtex_gene_medians, sym)
                for sym in high_id_paralog_syms
                for _ in ([0] if ref_gtex_tissue else [])
            }

            domains = (
                preloaded_pfam_domains
                if can_use_fast_domains
                else (fut_doms.result() if fut_doms is not None else [])
            )
            if fut_doms is not None or fut_ens_paralogs is not None:
                used_ensembl_rest = True
            flanked_seqs = fut_flanks.result() if fut_flanks is not None else None
            if fut_ccle_late is not None:
                ccle_profile = fut_ccle_late.result()
            if (
                flanked_seqs
                and chrom.replace("chr", "").upper() in GRCH38_NC_MAP
                and "NCBI E-utilities (GRCh38 RefSeq)" not in sources_used
            ):
                sources_used.append("NCBI E-utilities (GRCh38 RefSeq)")
            raw_paralogs = ot_paralogs if ot_paralogs is not None else (
                fut_ens_paralogs.result() if fut_ens_paralogs is not None else []
            )
            paralog_tpms: dict[str, Optional[float]] = {}
            for sym, p_fut in paralog_gtex_futs.items():
                p_meds = p_fut.result()
                if ref_gtex_tissue:
                    paralog_tpms[sym] = p_meds.get(ref_gtex_tissue)

            # Reconstruct exact 5'->3' CDS sequence from flanked_seqs (saving an extra Ensembl round-trip!)
            cds_seq = ""
            if flanked_seqs and all(flanked_seqs.get(iv) for iv in coding_intervals):
                cds_chunks: list[str] = []
                for ovs, ove in coding_intervals:
                    fseq = flanked_seqs[(ovs, ove)] or ""
                    lpad = min(flank_pad, ovs - 1)
                    clen = ove - ovs + 1
                    chunk = fseq[lpad : lpad + clen]
                    cds_chunks.append(chunk if strand >= 0 else reverse_complement(chunk))
                cds_seq = "".join(cds_chunks)
            elif fut_cds_fallback is not None:
                used_ensembl_rest = True
                cds_seq = fut_cds_fallback.result()
            else:
                used_ensembl_rest = True
                cds_seq = self.ensembl.fetch_cds_sequence(primary_summary.transcript_id)

        if used_ensembl_rest and "Ensembl REST API (GRCh38)" not in sources_used:
            sources_used.insert(0, "Ensembl REST API (GRCh38)")

        # If Open Targets paralogs were unavailable and fallback Ensembl paralogs found >=50% id hits:
        if ot_paralogs is None and raw_paralogs and ref_gtex_tissue:
            extra_syms = [
                str(p["gene_symbol"])
                for p in raw_paralogs[:3]
                if max(float(p["percent_identity"]), float(p["target_percent_identity"])) >= 50.0
            ][:2]
            for sym in extra_syms:
                p_meds = self.expression.fetch_gtex_gene_medians(sym)
                paralog_tpms[sym] = p_meds.get(ref_gtex_tissue)

        exons = score_transcript_exons(
            chrom=chrom,
            primary_tx=primary_raw,
            alt_coding_txs=alt_raws,
            domains=domains,
            cds_seq=cds_seq,
            flanked_seqs=flanked_seqs,
            flank_pad=flank_pad,
            gene_symbol=gene_symbol,
        )
        recommended_guides = select_top_gene_guides(exons, max_total=6)

        crispri_guides: list[GuideCandidate] = []
        if tss_pos and tss_iv and flanked_seqs and flanked_seqs.get(tss_iv):
            tss_win_start = max(1, tss_iv[0] - flank_pad)
            crispri_guides = design_crispri_tss_guides(
                tss_window_plus_seq=flanked_seqs[tss_iv] or "",
                window_genomic_start=tss_win_start,
                tss_pos=tss_pos,
                chrom=chrom,
                gene_strand=strand,
                gene_symbol=gene_symbol,
                max_guides=5,
            )

        if ccle_profile and ccle_profile.matched_ccle_sample:
            sources_used.append("cBioPortal CCLE (Broad DepMap Somatic Mutations & CNA)")
        if mygene_tss:
            sources_used.append("MyGene.info v3 (RefSeq GRCh38 TSS)")

        # 4. Essentiality & gnomAD constraint (only call direct gnomAD API if Open Targets lacked LOEUF)
        if essentiality.depmap_total_lines > 0 or ot_constraint or chem_tools:
            sources_used.append("Open Targets Platform v4 (DepMap Chronos, Tractability & Chemical Probes)")

        if ot_constraint:
            essentiality.gnomad_pli = ot_constraint.get("pli")
            essentiality.gnomad_loeuf = ot_constraint.get("loeuf")
            essentiality.gnomad_oe_lof = ot_constraint.get("oe_lof")

        if self.query_gnomad_direct and essentiality.gnomad_loeuf is None:
            g_const = self.gnomad.fetch_constraint(gene_symbol)
            if g_const:
                sources_used.append("gnomAD v4 GraphQL API")
                if g_const.get("pli") is not None:
                    essentiality.gnomad_pli = round(float(g_const["pli"]), 3)
                if g_const.get("loeuf") is not None:
                    essentiality.gnomad_loeuf = round(float(g_const["loeuf"]), 3)
                if g_const.get("oe_lof") is not None:
                    essentiality.gnomad_oe_lof = round(float(g_const["oe_lof"]), 3)

        # 5. Baseline tissue expression & HPA protein localization
        expr_profile = ExpressionProfile(
            requested_tissue=tissue,
            resolved_gtex_tissue=resolved_tid,
            requested_tissue_tpm=req_tpm,
            top_tissues=top_sorted,
            hpa_subcellular_locations=hpa_info["subcellular_locations"],
            hpa_tissue_ihc_level=hpa_info["ihc_level"],
            hpa_rna_specificity=hpa_info["rna_specificity"],
        )
        if target_medians:
            sources_used.append("GTEx Portal API v2 (v10 RNA-seq)")
        if expr_profile.hpa_subcellular_locations or expr_profile.hpa_rna_specificity:
            sources_used.append("Human Protein Atlas (HPA)")

        paralog_risks = classify_paralog_risks(raw_paralogs[:6], paralog_tpms)

        # 7. Modality verdict & trap warnings
        modality = evaluate_modality(
            gene_symbol=gene_symbol,
            primary_isoform=primary_summary,
            alt_isoforms=alt_summaries,
            essentiality=essentiality,
            expression=expr_profile,
            paralogs=paralog_risks,
            exons=exons,
            cell_line_genetics=ccle_profile,
        )

        # 8. Orthogonal chemical tools & wet-lab validation plan
        validation = build_validation_plan(
            gene_symbol=gene_symbol,
            primary_isoform=primary_summary,
            exons=exons,
            expression=expr_profile,
            chemical_tools=chem_tools,
        )

        now_iso = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        return PreflightReport(
            gene_symbol=gene_symbol,
            ensembl_id=ensembl_id,
            gene_description=raw_desc,
            chrom=f"chr{chrom.replace('chr', '')}",
            strand=strand,
            assembly="GRCh38",
            primary_transcript=primary_summary,
            alternative_isoforms=alt_summaries,
            essentiality=essentiality,
            expression=expr_profile,
            paralogs=paralog_risks,
            modality=modality,
            domains=domains,
            exons=exons,
            validation=validation,
            data_sources_used=sources_used,
            generated_at=now_iso,
            recommended_guides=recommended_guides,
            crispri_guides=crispri_guides,
            tss_genomic_pos=tss_pos,
            cell_line_genetics=ccle_profile,
        )
