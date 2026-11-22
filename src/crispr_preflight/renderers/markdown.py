"""Markdown renderer for Benchling, Notion, and GitHub lab notebooks."""

from __future__ import annotations

from crispr_preflight.models import PreflightReport


def render_markdown(report: PreflightReport) -> str:
    """Render a complete lab-notebook Markdown dossier from a PreflightReport."""
    tx = report.primary_transcript
    strand_label = "`+` forward" if report.strand >= 0 else "`-` reverse"
    mane_label = "MANE Select" if tx.is_mane_select else ("Ensembl Canonical" if tx.is_canonical else "Primary")

    nuc_lbl = report.nuclease_label or "SpCas9"
    pam_lbl = report.pam_label or "NGG"
    sp_len = report.spacer_length_nt or 20
    fmt_lbl = f"5'-{pam_lbl} + {sp_len}-nt" if report.nuclease_id == "cas12a" else f"{sp_len}-nt + {pam_lbl}"

    md: list[str] = [
        f"# CRISPR Preflight Dossier: `{report.gene_symbol}` (`{report.ensembl_id}`)",
        "",
        f"- **Organism:** *{report.species_scientific_name}* ({report.species_common_name})",
        f"- **CRISPR Nuclease:** **{nuc_lbl}** (`{fmt_lbl}`)",
        f"- **Description:** {report.gene_description or 'N/A'}",
        f"- **Genomic Locus ({report.assembly}):** `{report.chrom}` ({strand_label} strand)",
        f"- **Primary Transcript:** `{tx.transcript_id}` ({mane_label}, {tx.num_exons} exons, {tx.protein_length_aa} aa, +{len(report.alternative_isoforms)} alt coding isoforms)",
        f"- **Recommended Modality:** **{report.modality.primary_modality}** (`{report.modality.risk_level}` risk)",
        f"- **Verdict:** {report.modality.summary_verdict}",
        "",
        "## 1. Essentiality & Baseline Expression",
        "",
    ]

    ess = report.essentiality
    if ess.matched_cell_line:
        md.append(
            f"- **Cell Line Fitness ({ess.matched_cell_line}):** DepMap Chronos = `{ess.matched_cell_line_chronos:+.2f}`"
        )
    clg = report.cell_line_genetics
    if clg and clg.matched_ccle_sample:
        muts = "; ".join(clg.mutations) if clg.mutations else "No coding somatic mutations"
        cna = clg.cna_label or "Copy-number unannotated"
        md.append(
            f"- **Cell Line Endogenous Genetics (`{clg.matched_ccle_sample}`, CCLE):** {muts} | {cna}"
        )
    if ess.depmap_total_lines > 0 and ess.depmap_pan_essential_fraction is not None:
        md.append(
            f"- **Pan-Cancer DepMap Essentiality:** `{ess.depmap_essential_lines}/{ess.depmap_total_lines}` cell lines (`{int(ess.depmap_pan_essential_fraction * 100)}%`) with Chronos <= -1.0 (median `{ess.depmap_median_chronos:+.2f}`)"
        )
    if ess.gnomad_loeuf is not None or ess.gnomad_pli is not None:
        md.append(
            f"- **gnomAD Constraint:** LOEUF = `{ess.gnomad_loeuf}` | pLI = `{ess.gnomad_pli}`"
        )

    exp = report.expression
    if exp.resolved_gtex_tissue and exp.requested_tissue_tpm is not None:
        md.append(
            f"- **Target Tissue Expression (`{exp.resolved_gtex_tissue}`):** `{exp.requested_tissue_tpm} TPM`"
        )
    if exp.top_tissues:
        tops = ", ".join(f"`{t}` ({v} TPM)" for t, v in exp.top_tissues[:3])
        md.append(f"- **Top GTEx Tissues:** {tops}")
    if exp.hpa_subcellular_locations:
        md.append(f"- **HPA Subcellular Localization:** {', '.join(exp.hpa_subcellular_locations)}")

    md.extend(["", "## 2. Biological Failure-Mode Alerts", ""])
    if not report.modality.warnings:
        md.append("No major knockout biological traps flagged.")
    else:
        for w in report.modality.warnings:
            md.append(f"### [{w.severity}] {w.title}")
            md.append(f"- **Mechanism:** {w.explanation}")
            md.append(f"- **Bench Action:** {w.bench_action}")
            md.append("")

    md.extend(
        [
            "## 3. Coding Exon Targetability Map",
            "",
            f"| Exon | Score | Tier | AA Span | Length (Symm) | Isoform Retention | {nuc_lbl} {pam_lbl} | NMD Status | Domains | Target Region ({report.assembly}) |",
            "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- | :--- | :--- |",
        ]
    )
    for ex in report.exons:
        if not ex.is_coding:
            continue
        sym = f"{ex.exon_length_bp} bp (div3)" if ex.is_symmetric else f"{ex.exon_length_bp} bp"
        pam = str(ex.spcas9_pam_count) if ex.spcas9_pam_count is not None else "-"
        doms = "; ".join(d.split(" (")[0] for d in ex.overlapping_domains[:2]) or "-"
        md.append(
            f"| Exon {ex.exon_number} | **{ex.priority_score}** | `{ex.tier}` | {ex.aa_start}–{ex.aa_end} | {sym} | {int(ex.isoform_coverage_fraction * 100)}% | {pam} | `{ex.nmd_status}` | {doms} | `{ex.crispick_target_region}` |"
        )

    if report.domains:
        md.extend(["", "### Conserved Protein Domains (Deduplicated N → C)", ""])
        for idx, dm in enumerate(report.domains, start=1):
            md.append(f"- **D{idx}: {dm.description}** (`{dm.source_db}:{dm.accession}`, aa {dm.aa_start}–{dm.aa_end})")

    if report.recommended_guides:
        md.extend(
            [
                "",
                f"## 4. Candidate {nuc_lbl} Spacer Sequences (`{fmt_lbl}`) & Genotyping Primers",
                "",
                f"> **Off-Target Safety Check:** Verify 0–2 mismatch genome-wide specificity (PAM = `{pam_lbl}`) in [Cas-OFFinder](http://www.rgenome.net/cas-offinder/) before ordering oligos.",
                "",
                f"| Guide ID | Exon & Cut AA | {sp_len}-nt Spacer (`5'->3'`) | PAM | Cut Pos ({report.assembly}) | GC% | Score | Golden Gate Fwd / Rev Oligos | TIDE / ICE PCR Primers (Fwd / Rev) |",
                "| :--- | :--- | :--- | :---: | :--- | :---: | :---: | :--- | :--- |",
            ]
        )
        for g in report.recommended_guides:
            pr = g.genotyping_primers
            pr_col = (
                f"Fwd `{pr.fwd_seq}` ({pr.fwd_tm:.1f}°C) / Rev `{pr.rev_seq}` ({pr.rev_tm:.1f}°C) — {pr.amplicon_bp} bp (seq w/ {pr.recommended_seq_primer})"
                if pr
                else "-"
            )
            md.append(
                f"| **{g.guide_id}** | Exon {g.exon_number} (aa {g.aa_cut}) | `{g.spacer_20nt}` | `{g.pam}` | `{g.chrom}:{g.genomic_cut_pos}` ({g.genomic_strand}) | {g.gc_percent}% | {g.sequence_score}/100 | `5'-{g.cloning_fwd_oligo}-3'` / `5'-{g.cloning_rev_oligo}-3'` | {pr_col} |"
            )

    if report.crispri_guides:
        tss_str = f" (`{report.chrom}:{report.tss_genomic_pos}`)" if report.tss_genomic_pos else ""
        md.extend(
            [
                "",
                f"## 4b. Promoter CRISPRi (dCas9-KRAB) Guides (TSS -50 to +300 bp){tss_str}",
                "",
                "| Guide ID | TSS Offset | 20-nt Spacer (`5'->3'`) | PAM | Target Pos (GRCh38) | GC% | Score | U6 Fwd Oligo | U6 Rev Oligo |",
                "| :--- | :---: | :--- | :---: | :--- | :---: | :---: | :--- | :--- |",
            ]
        )
        for g in report.crispri_guides:
            off_s = f"{g.tss_offset_bp:+d} bp" if g.tss_offset_bp is not None else "-"
            md.append(
                f"| **{g.guide_id}** | `{off_s}` | `{g.spacer_20nt}` | `{g.pam}` | `{g.chrom}:{g.genomic_cut_pos}` ({g.genomic_strand}) | {g.gc_percent}% | {g.sequence_score}/100 | `5'-{g.cloning_fwd_oligo}-3'` | `5'-{g.cloning_rev_oligo}-3'` |"
            )

    if report.paralogs:
        md.extend(
            [
                "",
                "## 5. Human Paralog Compensation Check",
                "",
                "| Paralog | Ensembl ID | Query % ID | Target % ID | Tissue TPM | Buffering Risk |",
                "| :--- | :--- | :---: | :---: | :---: | :---: |",
            ]
        )
        for p in report.paralogs:
            tpm_s = f"{p.target_tissue_tpm}" if p.target_tissue_tpm is not None else "-"
            md.append(
                f"| **{p.gene_symbol}** | `{p.ensembl_id}` | {p.percent_identity:.1f}% | {p.target_percent_identity:.1f}% | {tpm_s} | `{p.buffering_risk}` |"
            )

    val = report.validation
    md.extend(
        [
            "",
            "## 6. Wet-Lab Validation & Orthogonal Controls",
            "",
            f"- **Recommended Protein Assay:** {val.recommended_protein_assay}",
            f"- **Antibody Epitope Guidance:** {val.epitope_guidance_summary}",
            f"- **RT-PCR Splicing Check:** {val.rt_qpcr_guidance}",
        ]
    )
    if val.orthogonal_chemical_tools:
        md.extend(["", "### Orthogonal Chemical Probes & Clinical Modulators", ""])
        for ct in val.orthogonal_chemical_tools:
            ctrl = " *(companion negative control available)*" if ct.has_negative_control else ""
            link = f" ([link]({ct.source_url}))" if ct.source_url else ""
            md.append(
                f"- **{ct.name}** (`{ct.tool_type}`, {ct.quality_or_stage}): {ct.mechanism_of_action}{ctrl}{link}"
            )

    md.extend(
        [
            "",
            "---",
            f"*Generated by **CRISPR Preflight** (Mulham Najajreh) at {report.generated_at}. Data sources: {', '.join(report.data_sources_used)}.*",
            "*This work was supported in part by complimentary access to Google Antigravity and Gemini 4 Argon to accelerate scientific workflows.*",
        ]
    )
    return "\n".join(md) + "\n"
