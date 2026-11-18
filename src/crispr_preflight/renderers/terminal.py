"""Rich zero-dependency ANSI terminal renderer for crispr-preflight reports."""

from __future__ import annotations

import os
import sys

from crispr_preflight.models import ExonAudit, PreflightReport


class _Style:
    def __init__(self, enabled: bool) -> None:
        self.on = enabled

    def _wrap(self, code: str, text: str) -> str:
        return f"\033[{code}m{text}\033[0m" if self.on else text

    def bold(self, t: str) -> str:
        return self._wrap("1", t)

    def dim(self, t: str) -> str:
        return self._wrap("2", t)

    def green(self, t: str) -> str:
        return self._wrap("1;32", t)

    def yellow(self, t: str) -> str:
        return self._wrap("1;33", t)

    def red(self, t: str) -> str:
        return self._wrap("1;31", t)

    def cyan(self, t: str) -> str:
        return self._wrap("1;36", t)


def render_terminal(
    report: PreflightReport,
    show_all_exons: bool = False,
    color: bool | None = None,
) -> str:
    """Format a PreflightReport for clear terminal inspection."""
    if color is None:
        color = sys.stdout.isatty() and ("NO_COLOR" not in os.environ)
    st = _Style(bool(color))
    lines: list[str] = []

    strand_str = "(+) forward" if report.strand >= 0 else "(-) reverse"
    tx = report.primary_transcript
    mane_badge = "MANE Select" if tx.is_mane_select else ("Canonical" if tx.is_canonical else "Primary")

    lines.append("=" * 92)
    lines.append(
        f" {st.bold('CRISPR-PREFLIGHT REPORT')}  |  {st.cyan(report.gene_symbol)} ({report.ensembl_id})  |  "
        f"{report.assembly} {report.chrom} {strand_str}"
    )
    if report.gene_description:
        lines.append(f" {st.dim(report.gene_description)}")
    lines.append(
        f" Primary Isoform: {st.bold(tx.transcript_id)} [{mane_badge}]  •  "
        f"{tx.num_exons} exon(s)  •  {tx.protein_length_aa} aa  •  "
        f"+{len(report.alternative_isoforms)} alt coding isoform(s)"
    )
    lines.append("=" * 92)

    # 1. Modality Verdict Banner
    mod = report.modality
    risk_badge = (
        st.red(f"[RISK: {mod.risk_level}]")
        if mod.risk_level == "HIGH"
        else (st.yellow(f"[RISK: {mod.risk_level}]") if mod.risk_level == "MODERATE" else st.green(f"[RISK: {mod.risk_level}]"))
    )
    lines.append(f"\n{st.bold('1. PERTURBATION STRATEGY VERDICT')}  {risk_badge}")
    lines.append(f"   Recommended Modality : {st.bold(mod.primary_modality)}")
    lines.append(f"   Summary              : {mod.summary_verdict}")

    # Essentiality & Expression quick metrics
    ess = report.essentiality
    ess_bits: list[str] = []
    if ess.matched_cell_line:
        ess_bits.append(f"{ess.matched_cell_line} Chronos = {ess.matched_cell_line_chronos:+.2f}")
    if ess.matched_lineage and ess.matched_lineage_median_chronos is not None:
        ess_bits.append(
            f"{ess.matched_lineage} median Chronos = {ess.matched_lineage_median_chronos:+.2f} (n={ess.matched_lineage_cell_count})"
        )
    if ess.depmap_total_lines > 0 and ess.depmap_pan_essential_fraction is not None:
        ess_bits.append(
            f"DepMap pan-essential: {ess.depmap_essential_lines}/{ess.depmap_total_lines} ({int(ess.depmap_pan_essential_fraction * 100)}%)"
        )
    if ess.gnomad_loeuf is not None:
        ess_bits.append(f"gnomAD LOEUF = {ess.gnomad_loeuf:.2f}")
    if ess.gnomad_pli is not None:
        ess_bits.append(f"pLI = {ess.gnomad_pli:.2f}")
    if ess_bits:
        lines.append(f"   Essentiality Metrics : {'  |  '.join(ess_bits)}")

    clg = report.cell_line_genetics
    if clg and clg.matched_ccle_sample:
        mut_s = "; ".join(clg.mutations) if clg.mutations else "Coding WT"
        cna_s = f"  |  {clg.cna_label}" if clg.cna_label else ""
        lines.append(f"   Cell-Line Genetics   : {clg.matched_ccle_sample} (CCLE): {mut_s}{cna_s}")

    exp = report.expression
    exp_bits: list[str] = []
    if exp.resolved_gtex_tissue and exp.requested_tissue_tpm is not None:
        exp_bits.append(f"{exp.resolved_gtex_tissue}: {exp.requested_tissue_tpm} TPM")
    if exp.top_tissues:
        top2 = ", ".join(f"{t} ({v} TPM)" for t, v in exp.top_tissues[:2])
        exp_bits.append(f"Top GTEx: {top2}")
    if exp.hpa_subcellular_locations:
        exp_bits.append(f"HPA Location: {', '.join(exp.hpa_subcellular_locations[:3])}")
    if exp_bits:
        lines.append(f"   Baseline Expression  : {'  |  '.join(exp_bits)}")

    # 2. Biological Trap Alerts
    lines.append(f"\n{st.bold('2. BIOLOGICAL TRAP ALERTS')}")
    if not mod.warnings:
        lines.append(f"   {st.green('✓ No major biological KO failure traps flagged.')}")
    else:
        for w in mod.warnings:
            sev_tag = (
                st.red(f"[{w.severity}]")
                if w.severity == "HIGH"
                else (st.yellow(f"[{w.severity}]") if w.severity == "MODERATE" else st.cyan(f"[{w.severity}]"))
            )
            lines.append(f"   {sev_tag} {st.bold(w.title)}")
            lines.append(f"       Why    : {w.explanation}")
            lines.append(f"       Action : {w.bench_action}")

    # 3. Exon-by-Exon Target Map
    coding_exons = [e for e in report.exons if e.is_coding]
    sorted_exons = sorted(coding_exons, key=lambda e: (-e.priority_score, e.exon_number))
    display_exons: list[ExonAudit] = coding_exons if show_all_exons else sorted_exons[:8]

    title_suffix = (
        f"(all {len(coding_exons)} coding exons in 5'->3' order)"
        if show_all_exons
        else f"(top {len(display_exons)} of {len(coding_exons)} coding exons ranked by safety; pass --all-exons for all)"
    )
    lines.append(f"\n{st.bold('3. TARGETABLE EXON AUDIT')} {st.dim(title_suffix)}")
    header = f"   {'Exon':<6} {'Score':<7} {'Tier':<13} {'AA Span':<11} {'Symm':<6} {'Isoforms':<9} {'NGG':<5} {'NMD Status':<20} {'CRISPick Region (GRCh38)'}"
    lines.append(st.bold(header))
    lines.append("   " + "-" * 104)

    for ex in display_exons:
        tier_str = (
            st.green(f"{ex.tier:<13}")
            if ex.tier == "RECOMMENDED"
            else (st.yellow(f"{ex.tier:<13}") if ex.tier == "CAUTION" else st.red(f"{ex.tier:<13}"))
        )
        aa_str = f"{ex.aa_start}-{ex.aa_end}" if ex.aa_start else "-"
        sym_str = "div3!" if ex.is_symmetric else "asym"
        iso_str = f"{int(ex.isoform_coverage_fraction * 100)}%"
        pam_str = str(ex.spcas9_pam_count) if ex.spcas9_pam_count is not None else "-"
        lines.append(
            f"   Ex {ex.exon_number:<3} {ex.priority_score:>3}/100 {tier_str} "
            f"{aa_str:<11} {sym_str:<6} {iso_str:<9} {pam_str:<5} "
            f"{ex.nmd_status:<20} {ex.crispick_target_region or '-'}"
        )
        # Show key domain / penalty notes concisely under top exons
        notes = [b for b in ex.score_breakdown[1:]]
        if notes:
            lines.append(f"          {st.dim('• ' + ' | '.join(notes[:2]))}")

    # 3b. Candidate SpCas9 sgRNA Spacer Sequences
    if report.recommended_guides:
        lines.append(
            f"\n{st.bold('4. CANDIDATE SpCas9 sgRNA SPACERS (20-nt + NGG)')} "
            f"{st.dim('(verify genome-wide off-targets in Cas-OFFinder / CRISPick before ordering)')}"
        )
        g_hdr = f"   {'Guide ID':<16} {'Exon/AA':<12} {'20-nt Spacer (5->3)':<22} {'PAM':<5} {'GC%':<5} {'Score':<7} {'U6 Fwd Oligo (BsmBI/BbsI)'}"
        lines.append(st.bold(g_hdr))
        lines.append("   " + "-" * 104)
        for g in report.recommended_guides:
            ex_aa = f"E{g.exon_number} aa{g.aa_cut}"
            lines.append(
                f"   {g.guide_id:<16} {ex_aa:<12} {st.cyan(g.spacer_20nt):<22} {g.pam:<5} "
                f"{str(g.gc_percent) + '%':<5} {str(g.sequence_score) + '/100':<7} 5'-{g.cloning_fwd_oligo}-3'"
            )
            if g.genotyping_primers:
                pr = g.genotyping_primers
                lines.append(
                    f"      {st.dim(f'TIDE/ICE PCR ({pr.amplicon_bp} bp, seq w/ {pr.recommended_seq_primer}): Fwd {pr.fwd_seq} ({pr.fwd_tm:.1f}C) | Rev {pr.rev_seq} ({pr.rev_tm:.1f}C)')}"
                )

    if report.crispri_guides:
        tss_info = f" (TSS {report.chrom}:{report.tss_genomic_pos})" if report.tss_genomic_pos else ""
        lines.append(f"\n{st.bold('4b. PROMOTER CRISPRi (dCas9-KRAB) GUIDES (TSS -50..+300 bp)')}{st.dim(tss_info)}")
        for g in report.crispri_guides[:4]:
            off_s = f"{g.tss_offset_bp:+d}bp" if g.tss_offset_bp is not None else "-"
            lines.append(
                f"   {g.guide_id:<18} {off_s:<9} {st.cyan(g.spacer_20nt)} {g.pam}  "
                f"{g.gc_percent}% GC ({g.sequence_score}/100)  5'-{g.cloning_fwd_oligo}-3'"
            )

    # 5. Paralog Buffering Table
    if report.paralogs:
        lines.append(f"\n{st.bold('5. HUMAN PARALOG BUFFERING CHECK')}")
        for p in report.paralogs[:4]:
            r_tag = (
                st.red(f"[{p.buffering_risk}]")
                if p.buffering_risk == "HIGH"
                else (st.yellow(f"[{p.buffering_risk}]") if p.buffering_risk == "MODERATE" else st.dim(f"[{p.buffering_risk}]"))
            )
            tpm_note = f", target tissue = {p.target_tissue_tpm} TPM" if p.target_tissue_tpm is not None else ""
            lines.append(
                f"   {r_tag} {st.bold(p.gene_symbol)} ({p.ensembl_id}) — "
                f"{p.percent_identity:.1f}% query id / {p.target_percent_identity:.1f}% target id{tpm_note}"
            )

    # 6. Bench Validation & Orthogonal Chemical Controls
    val = report.validation
    lines.append(f"\n{st.bold('6. WET-LAB VALIDATION & ORTHOGONAL CONTROLS')}")
    lines.append(f"   Protein Assay     : {val.recommended_protein_assay}")
    lines.append(f"   Antibody Epitopes : {val.epitope_guidance_summary}")
    avoid_wins = [w for w in val.epitope_windows if w.status == "AVOID"]
    for aw in avoid_wins[:2]:
        lines.append(f"       {st.yellow('Avoid aa ' + str(aw.aa_start) + '-' + str(aw.aa_end))}: {aw.reason}")
    lines.append(f"   RT-PCR Check      : {val.rt_qpcr_guidance}")

    if val.orthogonal_chemical_tools:
        lines.append("   Orthogonal Small-Molecule Probes / Modulators:")
        for ct in val.orthogonal_chemical_tools[:4]:
            ctrl_note = " [+neg control available]" if ct.has_negative_control else ""
            lines.append(
                f"       • {st.bold(ct.name)} ({ct.tool_type}, {ct.quality_or_stage}) — {ct.mechanism_of_action}{ctrl_note}"
            )
    else:
        lines.append("   Orthogonal Probes : None curated in Open Targets (rely on cDNA rescue / multiple independent sgRNAs).")

    lines.append("\n" + "-" * 92)
    lines.append(
        f" {st.dim('Data sources: ' + ', '.join(report.data_sources_used) + '. Review respective database terms for reuse.')}"
    )
    return "\n".join(lines) + "\n"
