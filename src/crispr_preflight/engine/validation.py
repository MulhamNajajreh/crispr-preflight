"""Wet-lab validation planner: antibody safe-epitope calculator, RT-PCR checks, and chemical probes."""

from __future__ import annotations

from crispr_preflight.models import (
    ChemicalTool,
    EpitopeWindow,
    ExonAudit,
    ExpressionProfile,
    IsoformSummary,
    ValidationPlan,
)


def build_validation_plan(
    gene_symbol: str,
    primary_isoform: IsoformSummary,
    exons: list[ExonAudit],
    expression: ExpressionProfile,
    chemical_tools: list[ChemicalTool],
) -> ValidationPlan:
    """Compute safe vs. misleading antibody epitope zones and recommend orthogonal bench controls."""
    prot_len = primary_isoform.protein_length_aa
    coding_exons = [e for e in exons if e.is_coding and e.aa_start and e.aa_end]

    # 1. Subcellular localization -> best protein readout assay
    locs_lower = " ".join(expression.hpa_subcellular_locations).lower()
    loc_label = ", ".join(expression.hpa_subcellular_locations) or "Unannotated localization"

    if "plasma membrane" in locs_lower or "cell junctions" in locs_lower:
        assay_rec = (
            f"Flow Cytometry (live surface staining) + Western Blot [{loc_label}] — "
            "surface FACS allows rapid single-cell KO quantification before clonal sorting."
        )
    elif "secreted" in locs_lower or "vesicles" in locs_lower:
        assay_rec = (
            f"Conditioned-Media ELISA / Western Blot [{loc_label}] — "
            "account for paracrine rescue from unedited neighbor cells in pooled cultures!"
        )
    elif "nucleoplasm" in locs_lower or "nuclear" in locs_lower or "chromatin" in locs_lower:
        assay_rec = (
            f"Western Blot (nuclear/RIPA whole-cell lysate) or Immunofluorescence [{loc_label}]."
        )
    else:
        assay_rec = f"Western Blot (whole-cell lysate) + Mass Spec / IF confirmation [{loc_label}]."

    # 2. Antibody Epitope Safe vs. Unsafe Amino-Acid Windows
    # Unsafe regions:
    #  (a) Pre-reinitiation N-terminal peptide (1 .. reinit_met) where truncated protein drops the N-term epitope!
    #  (b) Symmetric skippable internal exons (if your sgRNA cuts Exon k and the antibody also binds Exon k,
    #      in-frame exon skipping gives a 100% blank Western blot on a still-expressed protein!)
    #  (c) Alternative-isoform skipped exons.
    avoid_windows: list[EpitopeWindow] = []
    max_reinit = max(
        (e.reinitiation_met_codon for e in coding_exons if e.reinitiation_met_codon),
        default=0,
    )
    if max_reinit > 1:
        avoid_windows.append(
            EpitopeWindow(
                aa_start=1,
                aa_end=max_reinit,
                status="AVOID",
                reason=(
                    f"Upstream of internal in-frame Met{max_reinit}: N-terminally truncated "
                    "re-initiation products will look 'blank' to N-terminal antibodies."
                ),
            )
        )

    for ex in coding_exons:
        assert ex.aa_start is not None and ex.aa_end is not None
        if ex.is_symmetric:
            avoid_windows.append(
                EpitopeWindow(
                    aa_start=ex.aa_start,
                    aa_end=ex.aa_end,
                    status="AVOID",
                    reason=(
                        f"Encoded by symmetric Exon {ex.exon_number} ({ex.exon_length_bp} bp): if cutting Exon {ex.exon_number} "
                        "causes in-frame exon skipping, an antibody binding here loses its epitope while internal-deletion protein persists."
                    ),
                )
            )
        elif not ex.is_constitutive:
            avoid_windows.append(
                EpitopeWindow(
                    aa_start=ex.aa_start,
                    aa_end=ex.aa_end,
                    status="AVOID",
                    reason=(
                        f"Encoded by non-constitutive Exon {ex.exon_number} (absent in {', '.join(ex.skipped_in_isoforms[:2])})."
                    ),
                )
            )

    # Compute complement constitutive SAFE windows downstream of N-term re-initiation
    unsafe_mask = [False] * (prot_len + 1)
    for w in avoid_windows:
        for aa in range(max(1, w.aa_start), min(prot_len, w.aa_end) + 1):
            unsafe_mask[aa] = True

    safe_windows: list[EpitopeWindow] = []
    run_start = None
    for aa in range(1, prot_len + 1):
        if not unsafe_mask[aa]:
            if run_start is None:
                run_start = aa
        else:
            if run_start is not None and (aa - run_start) >= 25:
                safe_windows.append(
                    EpitopeWindow(
                        aa_start=run_start,
                        aa_end=aa - 1,
                        status="SAFE",
                        reason="Constitutive region outside N-terminal re-initiation and symmetric skippable exons.",
                    )
                )
            run_start = None
    if run_start is not None and (prot_len - run_start + 1) >= 25:
        safe_windows.append(
            EpitopeWindow(
                aa_start=run_start,
                aa_end=prot_len,
                status="SAFE",
                reason="Constitutive C-terminal/core region downstream of re-initiation and skippable exons.",
            )
        )

    all_windows = safe_windows[:3] + avoid_windows[:4]
    if safe_windows:
        best_safe = max(safe_windows, key=lambda w: w.aa_end - w.aa_start)
        epitope_summary = (
            f"Choose a Western antibody recognizing residues aa {best_safe.aa_start}–{best_safe.aa_end} "
            "(never rely solely on an antibody whose immunogen overlaps your sgRNA-targeted exon)."
        )
    else:
        epitope_summary = (
            "Use polyclonal or two distinct monoclonal antibodies (N- and C-terminal) to rule out truncated escapers."
        )

    top_exons = [e for e in coding_exons if e.tier == "RECOMMENDED"]
    if top_exons:
        e0 = top_exons[0]
        rt_guidance = (
            f"Design RT-PCR primers in exons flanking Exon {e0.exon_number} (not inside the cut site!) "
            f"to detect both NMD transcript knockdown and any unexpected exon-skipping band shift."
        )
    else:
        rt_guidance = (
            "Design RT-PCR primers spanning exons upstream and downstream of the edited exon to check for cryptic splicing."
        )

    return ValidationPlan(
        recommended_protein_assay=assay_rec,
        epitope_guidance_summary=epitope_summary,
        epitope_windows=all_windows,
        rt_qpcr_guidance=rt_guidance,
        orthogonal_chemical_tools=chemical_tools,
    )
