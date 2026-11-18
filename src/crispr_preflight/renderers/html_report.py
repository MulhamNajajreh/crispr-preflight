"""Standalone, zero-dependency interactive HTML + SVG Exon/Domain/NMD report renderer."""

from __future__ import annotations

import html
import json

from crispr_preflight.models import PreflightReport


def _build_architecture_svg(report: PreflightReport) -> str:
    """Generate a proportional 5'->3' protein/exon/domain/NMD SVG track diagram."""
    prot_len = max(1, report.primary_transcript.protein_length_aa)
    width = 920
    pad_x = 45
    usable_w = width - 2 * pad_x

    def aa_to_x(aa: float) -> float:
        clamped = max(1.0, min(float(prot_len), aa))
        return pad_x + ((clamped - 1.0) / max(1.0, float(prot_len - 1))) * usable_w

    svg: list[str] = [
        f'<svg viewBox="0 0 {width} 225" width="100%" role="img" aria-label="Exon, Domain, and NMD Track Map" style="font-family: system-ui, -apple-system, sans-serif; background:#ffffff; border:1px solid #e2e8f0; border-radius:10px;">',
        f'<text x="{pad_x}" y="24" font-size="12" font-weight="600" fill="#334155">1. Coding Exons (5\' → 3\' colored by CRISPR KO Safety Tier — click any exon row below for score math)</text>',
    ]

    tier_colors = {
        "RECOMMENDED": ("#10b981", "#047857"),
        "CAUTION": ("#f59e0b", "#b45309"),
        "AVOID": ("#ef4444", "#b91c1c"),
    }

    # Track 1: Coding Exons proportional to AA span
    svg.append(f'<line x1="{pad_x}" y1="54" x2="{width - pad_x}" y2="54" stroke="#cbd5e1" stroke-width="3"/>')
    for ex in report.exons:
        if not ex.is_coding or not ex.aa_start or not ex.aa_end:
            continue
        x1 = aa_to_x(ex.aa_start)
        x2 = max(x1 + 5.0, aa_to_x(ex.aa_end))
        fill, stroke = tier_colors.get(ex.tier, ("#94a3b8", "#475569"))
        sym_tag = " (symmetric div3)" if ex.is_symmetric else ""
        tip = html.escape(
            f"Exon {ex.exon_number}: Score {ex.priority_score}/100 [{ex.tier}] | aa {ex.aa_start}-{ex.aa_end}{sym_tag} | NMD: {ex.nmd_status}"
        )
        svg.append(
            f'<g><title>{tip}</title>'
            f'<rect x="{x1:.1f}" y="38" width="{(x2 - x1):.1f}" height="32" rx="4" fill="{fill}" fill-opacity="0.88" stroke="{stroke}" stroke-width="1.2"/>'
        )
        if (x2 - x1) >= 24:
            svg.append(
                f'<text x="{(x1 + x2) / 2:.1f}" y="58" text-anchor="middle" font-size="10" font-weight="700" fill="#ffffff">E{ex.exon_number}</text>'
            )
        svg.append("</g>")

    # Track 2: Conserved Pfam / Structural Domains (multi-lane non-overlapping layout + full Domain Key below)
    lane_ends: list[float] = []
    dom_elements: list[str] = []
    for idx, dom in enumerate(report.domains, start=1):
        dx1 = aa_to_x(dom.aa_start)
        dx2 = max(dx1 + 18.0, aa_to_x(dom.aa_end))
        lane_idx = -1
        for li, end_x in enumerate(lane_ends):
            if dx1 >= end_x + 4.0:
                lane_idx = li
                lane_ends[li] = dx2
                break
        if lane_idx == -1:
            lane_idx = len(lane_ends)
            lane_ends.append(dx2)

        dy = 108 + lane_idx * 28
        dw = dx2 - dx1
        dtip = html.escape(f"D{idx}: {dom.description} [{dom.source_db}:{dom.accession}] (aa {dom.aa_start}-{dom.aa_end})")
        full_lbl = f"D{idx}: {dom.description}"
        bar_lbl = full_lbl if dw >= (len(full_lbl) * 6.8 + 14.0) else (f"D{idx}" if dw >= 22.0 else "")
        dom_elements.append(
            f'<g><title>{dtip}</title>'
            f'<rect x="{dx1:.1f}" y="{dy}" width="{dw:.1f}" height="24" rx="4" fill="#4f46e5" fill-opacity="0.9" stroke="#312e81" stroke-width="1"/>'
        )
        if bar_lbl:
            dom_elements.append(
                f'<text x="{(dx1 + dx2) / 2:.1f}" y="{dy + 16}" text-anchor="middle" font-size="11" font-weight="700" fill="#ffffff">{html.escape(bar_lbl)}</text>'
            )
        dom_elements.append("</g>")

    num_lanes = max(1, len(lane_ends))
    dom_track_h = num_lanes * 28 - 4
    t3_label_y = 108 + dom_track_h + 24
    t3_rect_y = t3_label_y + 10
    axis_y = t3_rect_y + 38
    svg_h = axis_y + 15
    svg[0] = (
        f'<svg viewBox="0 0 {width} {svg_h}" width="100%" role="img" aria-label="Exon, Domain, and NMD Track Map" '
        f'style="font-family: system-ui, -apple-system, sans-serif; background:#ffffff; border:1px solid #e2e8f0; border-radius:10px;">'
    )

    svg.append(
        f'<text x="{pad_x}" y="98" font-size="12" font-weight="600" fill="#334155">2. Conserved Protein Domains (cutting here neutralizes -3 bp in-frame escapers)</text>'
    )
    svg.append(f'<rect x="{pad_x}" y="108" width="{usable_w}" height="{dom_track_h}" rx="4" fill="#f1f5f9"/>')
    if dom_elements:
        svg.extend(dom_elements)
    else:
        svg.append(
            f'<text x="{width / 2:.1f}" y="124" text-anchor="middle" font-size="11" fill="#64748b">No Pfam/Smart domains annotated on this isoform</text>'
        )

    # Track 3: Antibody Epitope Safe vs Avoid Windows
    svg.append(
        f'<text x="{pad_x}" y="{t3_label_y}" font-size="12" font-weight="600" fill="#334155">3. Western Antibody Epitope Map (Teal = Safe Constitutive Epitope | Coral = Re-initiation / Skippable Exon Risk)</text>'
    )
    svg.append(f'<rect x="{pad_x}" y="{t3_rect_y}" width="{usable_w}" height="18" rx="4" fill="#e2e8f0"/>')
    for win in report.validation.epitope_windows:
        wx1 = aa_to_x(win.aa_start)
        wx2 = max(wx1 + 4.0, aa_to_x(win.aa_end))
        wfill = "#0d9488" if win.status == "SAFE" else "#f43f5e"
        wtip = html.escape(f"{win.status} (aa {win.aa_start}-{win.aa_end}): {win.reason}")
        svg.append(
            f'<g><title>{wtip}</title><rect x="{wx1:.1f}" y="{t3_rect_y}" width="{(wx2 - wx1):.1f}" height="18" rx="3" fill="{wfill}" fill-opacity="0.82"/></g>'
        )

    # Axis ticks
    svg.append(f'<text x="{pad_x}" y="{axis_y}" font-size="11" fill="#64748b">N-term (aa 1)</text>')
    svg.append(
        f'<text x="{width / 2:.1f}" y="{axis_y}" text-anchor="middle" font-size="11" fill="#64748b">aa {prot_len // 2}</text>'
    )
    svg.append(
        f'<text x="{width - pad_x}" y="{axis_y}" text-anchor="end" font-size="11" fill="#64748b">C-term (aa {prot_len})</text>'
    )
    svg.append("</svg>")
    return "\n".join(svg)


def render_html(report: PreflightReport) -> str:
    """Render a self-contained HTML lab-notebook report with interactive exon score breakdown & CRISPick copy."""
    esc = html.escape
    tx = report.primary_transcript
    mod = report.modality
    ess = report.essentiality
    exp = report.expression
    val = report.validation

    risk_bg = {"HIGH": "#fef2f2", "MODERATE": "#fffbeb", "LOW": "#ecfdf5"}.get(mod.risk_level, "#f8fafc")
    risk_fg = {"HIGH": "#b91c1c", "MODERATE": "#b45309", "LOW": "#047857"}.get(mod.risk_level, "#334155")
    risk_border = {"HIGH": "#fecaca", "MODERATE": "#fde68a", "LOW": "#a7f3d0"}.get(mod.risk_level, "#cbd5e1")

    svg_diagram = _build_architecture_svg(report)

    warnings_html = ""
    if not mod.warnings:
        warnings_html = '<div class="card ok-card">✓ No high-severity knockout failure traps detected.</div>'
    else:
        cards: list[str] = []
        for w in mod.warnings:
            cls = "warn-high" if w.severity == "HIGH" else ("warn-mod" if w.severity == "MODERATE" else "warn-info")
            cards.append(
                f'<div class="trap-card {cls}">'
                f'<div class="trap-hdr"><span class="badge">{esc(w.severity)}</span> <strong>{esc(w.title)}</strong></div>'
                f'<p><strong>Biology:</strong> {esc(w.explanation)}</p>'
                f'<p><strong>Bench Action:</strong> {esc(w.bench_action)}</p>'
                f"</div>"
            )
        warnings_html = "\n".join(cards)

    exon_rows: list[str] = []
    for ex in report.exons:
        if not ex.is_coding:
            continue
        t_cls = (
            "tier-rec"
            if ex.tier == "RECOMMENDED"
            else ("tier-caut" if ex.tier == "CAUTION" else "tier-avoid")
        )
        dom_txt = ", ".join(d.split(" (")[0] for d in ex.overlapping_domains[:2]) or "—"
        sym_txt = f"{ex.exon_length_bp} bp (div3)" if ex.is_symmetric else f"{ex.exon_length_bp} bp"
        pam_txt = str(ex.spcas9_pam_count) if ex.spcas9_pam_count is not None else "—"
        reg = ex.crispick_target_region or ""
        bd_items = "".join(f"<li>{esc(item)}</li>" for item in ex.score_breakdown)
        exon_rows.append(
            f"<tr>"
            f"<td><strong>Exon {ex.exon_number}</strong></td>"
            f"<td><strong>{ex.priority_score}</strong>/100</td>"
            f'<td><span class="pill {t_cls}">{esc(ex.tier)}</span></td>'
            f"<td>aa {ex.aa_start}–{ex.aa_end}</td>"
            f"<td>{esc(sym_txt)}</td>"
            f"<td>{int(ex.isoform_coverage_fraction * 100)}%</td>"
            f"<td>{esc(pam_txt)}</td>"
            f"<td><code>{esc(ex.nmd_status)}</code></td>"
            f"<td>{esc(dom_txt)}</td>"
            f'<td><code>{esc(reg)}</code> <button class="copy-btn" data-copy="{esc(reg)}" onclick="copyCoord(this)">Copy</button></td>'
            f"</tr>"
            f'<tr class="breakdown-row"><td colspan="10"><ul class="bd-list">{bd_items}</ul></td></tr>'
        )

    domain_key_html = ""
    if report.domains:
        chips = []
        for idx, dm in enumerate(report.domains, start=1):
            chips.append(
                f'<span style="display:inline-block;background:#ffffff;border:1px solid #a5b4fc;border-radius:6px;padding:3px 8px;margin:3px 4px 3px 0;font-size:12px;">'
                f'<strong style="color:#4338ca;">D{idx}</strong>: <strong>{esc(dm.description)}</strong> '
                f'<span style="color:#475569;">({esc(dm.source_db)}:{esc(dm.accession)}, aa {dm.aa_start}–{dm.aa_end})</span></span>'
            )
        domain_key_html = (
            '<div style="background:#eef2ff;border:1px solid #c7d2fe;border-radius:8px;padding:8px 12px;margin-top:8px;">'
            '<div style="font-size:11.5px;font-weight:700;color:#312e81;margin-bottom:4px;">Track 2 Conserved Domain Key (complete names):</div>'
            f'{"".join(chips)}</div>'
        )

    guides_section = ""
    if report.recommended_guides:
        g_rows = []
        for g in report.recommended_guides:
            dom_tag = f" ({esc(g.overlaps_domain)})" if g.overlaps_domain else ""
            pr = g.genotyping_primers
            if pr:
                primer_html = (
                    f"<div>Fwd ({pr.fwd_tm:.1f}&deg;C): <code>5'-{esc(pr.fwd_seq)}-3'</code></div>"
                    f"<div>Rev ({pr.rev_tm:.1f}&deg;C): <code>5'-{esc(pr.rev_seq)}-3'</code></div>"
                    f'<div style="font-size:11px;color:#475569;">{pr.amplicon_bp} bp &bull; Seq w/ {esc(pr.recommended_seq_primer)}</div>'
                )
            else:
                primer_html = "—"
            g_rows.append(
                f"<tr><td><strong>{esc(g.guide_id)}</strong></td>"
                f"<td>Exon {g.exon_number} (aa {g.aa_cut}){dom_tag}</td>"
                f'<td><code>{esc(g.spacer_20nt)}</code> <strong>{esc(g.pam)}</strong> <button class="copy-btn" data-copy="{esc(g.spacer_20nt)}" onclick="copyCoord(this)">Copy</button></td>'
                f"<td><code>{esc(g.chrom)}:{g.genomic_cut_pos}</code> ({esc(g.genomic_strand)} / {esc(g.transcript_orientation)})</td>"
                f"<td>{g.gc_percent}% GC ({g.sequence_score}/100)</td>"
                f"<td><code>5'-{esc(g.cloning_fwd_oligo)}-3'</code><br/><code>5'-{esc(g.cloning_rev_oligo)}-3'</code></td>"
                f"<td>{primer_html}</td></tr>"
            )
        guides_section = (
            "<h2>4. Candidate SpCas9 sgRNA Spacer Sequences (20-nt + NGG), U6 Oligos &amp; TIDE/ICE Primers</h2>"
            '<p style="font-size:13px;color:#1e3a8a;background:#eff6ff;border:1px solid #bfdbfe;padding:9px 12px;border-radius:7px;">'
            "<strong>Off-target verification step:</strong> These 20-nt spacers are filtered for U6 Pol III compatibility (no <code>TTTT</code>), GC balance, Graf PAM-proximal efficiency, and BsmBI/BbsI cloning overhangs. "
            "Verify zero 0–2 mismatch genomic off-targets in <a href=\"http://www.rgenome.net/cas-offinder/\" target=\"_blank\" rel=\"noopener\">Cas-OFFinder</a> or "
            "<a href=\"https://portals.broadinstitute.org/gppx/crispick/public\" target=\"_blank\" rel=\"noopener\">Broad CRISPick</a> before ordering oligos.</p>"
            f"<table><thead><tr><th>Guide ID</th><th>Exon &amp; AA</th><th>20-nt Spacer + PAM</th><th>Cut Pos ({esc(report.assembly)})</th><th>GC% &amp; Score</th><th>U6 Golden Gate Oligos</th><th>TIDE / Synthego ICE PCR Primers</th></tr></thead>"
            f"<tbody>{''.join(g_rows)}</tbody></table>"
        )

    crispri_section = ""
    if report.crispri_guides:
        ci_rows = []
        for g in report.crispri_guides:
            off_s = f"{g.tss_offset_bp:+d} bp" if g.tss_offset_bp is not None else "—"
            ci_rows.append(
                f"<tr><td><strong>{esc(g.guide_id)}</strong></td>"
                f"<td><code>{esc(off_s)}</code> from TSS</td>"
                f'<td><code>{esc(g.spacer_20nt)}</code> <strong>{esc(g.pam)}</strong> <button class="copy-btn" data-copy="{esc(g.spacer_20nt)}" onclick="copyCoord(this)">Copy</button></td>'
                f"<td><code>{esc(g.chrom)}:{g.genomic_cut_pos}</code> ({esc(g.genomic_strand)})</td>"
                f"<td>{g.gc_percent}% GC ({g.sequence_score}/100)</td>"
                f"<td><code>5'-{esc(g.cloning_fwd_oligo)}-3'</code><br/><code>5'-{esc(g.cloning_rev_oligo)}-3'</code></td></tr>"
            )
        tss_lbl = f" ({esc(report.assembly)} <code>{esc(report.chrom)}:{report.tss_genomic_pos}</code>)" if report.tss_genomic_pos else ""
        crispri_section = (
            f"<h2>4b. Promoter CRISPRi (dCas9-KRAB) Candidate Guides (TSS -50 to +300 bp){tss_lbl}</h2>"
            f"<table><thead><tr><th>Guide ID</th><th>TSS Offset</th><th>20-nt Spacer + PAM</th><th>Target Pos ({esc(report.assembly)})</th><th>GC% &amp; Score</th><th>U6 Golden Gate Oligos</th></tr></thead>"
            f"<tbody>{''.join(ci_rows)}</tbody></table>"
        )

    paralog_section = ""
    if report.paralogs:
        p_rows = []
        for p in report.paralogs:
            tpm_v = f"{p.target_tissue_tpm} TPM" if p.target_tissue_tpm is not None else "—"
            p_rows.append(
                f"<tr><td><strong>{esc(p.gene_symbol)}</strong> (<code>{esc(p.ensembl_id)}</code>)</td>"
                f"<td>{p.percent_identity:.1f}%</td><td>{p.target_percent_identity:.1f}%</td>"
                f"<td>{esc(tpm_v)}</td><td><strong>{esc(p.buffering_risk)}</strong></td></tr>"
            )
        paralog_section = (
            f"<h2>5. {esc(report.species_common_name)} Paralog Buffering Check</h2>"
            "<table><thead><tr><th>Paralog</th><th>Query % ID</th><th>Target % ID</th><th>Tissue Expression</th><th>Buffering Risk</th></tr></thead>"
            f"<tbody>{''.join(p_rows)}</tbody></table>"
        )

    chem_items = ""
    if val.orthogonal_chemical_tools:
        lis = []
        for ct in val.orthogonal_chemical_tools:
            ctrl = ' <span class="pill tier-rec">+Neg Control</span>' if ct.has_negative_control else ""
            link = f' — <a href="{esc(ct.source_url)}" target="_blank" rel="noopener">Open Resource</a>' if ct.source_url else ""
            lis.append(
                f"<li><strong>{esc(ct.name)}</strong> ({esc(ct.tool_type)}, {esc(ct.quality_or_stage)}): "
                f"{esc(ct.mechanism_of_action)}{ctrl}{link}</li>"
            )
        chem_items = f"<ul>{''.join(lis)}</ul>"
    else:
        chem_items = "<p>No curated selective chemical probes listed in Open Targets; rely on cDNA rescue and >=2 independent sgRNAs.</p>"

    chronos_metric = "—"
    if ess.matched_cell_line and ess.matched_cell_line_chronos is not None:
        chronos_metric = f"{ess.matched_cell_line}: {ess.matched_cell_line_chronos:+.2f}"
    elif ess.depmap_median_chronos is not None:
        chronos_metric = f"Median: {ess.depmap_median_chronos:+.2f} ({int((ess.depmap_pan_essential_fraction or 0) * 100)}% pan-ess)"

    loeuf_metric = f"{ess.gnomad_loeuf:.2f}" if ess.gnomad_loeuf is not None else "—"
    tpm_metric = (
        f"{exp.requested_tissue_tpm} TPM ({exp.resolved_gtex_tissue})"
        if exp.resolved_gtex_tissue and exp.requested_tissue_tpm is not None
        else (f"{exp.top_tissues[0][1]} TPM ({exp.top_tissues[0][0]})" if exp.top_tissues else "—")
    )

    ccle_sub = ""
    clg = report.cell_line_genetics
    if clg and clg.matched_ccle_sample:
        mut_txt = ", ".join(clg.mutations[:2]) if clg.mutations else "Coding WT"
        cna_txt = f" • {clg.cna_label}" if clg.cna_label else ""
        ccle_sub = f'<div style="font-size:11.5px;color:#475569;margin-top:3px;">CCLE ({esc(clg.matched_ccle_sample)}): {esc(mut_txt)}{esc(cna_txt)}</div>'

    json_payload = json.dumps(report.to_dict()).replace("</", "<\\/")

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>CRISPR Preflight — {esc(report.gene_symbol)}</title>
<style>
  body {{ font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; color:#0f172a; background:#f8fafc; margin:0; padding:24px; line-height:1.5; }}
  .container {{ max-width: 1120px; margin: 0 auto; background:#ffffff; border-radius:14px; padding:28px 32px; box-shadow:0 4px 20px rgba(15,23,42,0.06); }}
  h1 {{ margin:0 0 6px 0; font-size:24px; }}
  h2 {{ margin-top:28px; font-size:18px; border-bottom:1px solid #e2e8f0; padding-bottom:6px; }}
  .sub {{ color:#475569; font-size:14px; margin-bottom:18px; }}
  .verdict {{ background:{risk_bg}; border:1px solid {risk_border}; color:{risk_fg}; border-radius:10px; padding:16px 20px; margin-bottom:20px; }}
  .kpi-grid {{ display:grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap:12px; margin-bottom:22px; }}
  .kpi {{ background:#f8fafc; border:1px solid #e2e8f0; border-radius:8px; padding:12px 14px; }}
  .kpi .lbl {{ font-size:12px; color:#64748b; text-transform:uppercase; letter-spacing:0.04em; }}
  .kpi .val {{ font-size:16px; font-weight:700; margin-top:4px; color:#0f172a; }}
  .trap-card {{ border-radius:8px; padding:12px 16px; margin-bottom:10px; border-left:4px solid #94a3b8; background:#f8fafc; }}
  .trap-card p {{ margin:4px 0; font-size:13.5px; }}
  .warn-high {{ border-left-color:#ef4444; background:#fef2f2; }}
  .warn-mod {{ border-left-color:#f59e0b; background:#fffbeb; }}
  .warn-info {{ border-left-color:#3b82f6; background:#eff6ff; }}
  .badge {{ font-size:11px; font-weight:700; padding:2px 6px; border-radius:4px; background:rgba(0,0,0,0.08); }}
  table {{ width:100%; border-collapse:collapse; font-size:13px; margin-top:10px; }}
  th, td {{ text-align:left; padding:9px 10px; border-bottom:1px solid #e2e8f0; }}
  th {{ background:#f1f5f9; font-weight:600; color:#334155; }}
  .breakdown-row td {{ background:#fcfcfd; padding:4px 18px 10px 28px; color:#475569; font-size:12px; }}
  .bd-list {{ margin:0; padding-left:16px; }}
  .pill {{ display:inline-block; padding:2px 8px; border-radius:999px; font-size:11px; font-weight:700; }}
  .tier-rec {{ background:#d1fae5; color:#065f46; }}
  .tier-caut {{ background:#fef3c7; color:#92400e; }}
  .tier-avoid {{ background:#fee2e2; color:#991b1b; }}
  .copy-btn {{ border:1px solid #cbd5e1; background:#ffffff; border-radius:5px; padding:2px 7px; font-size:11px; cursor:pointer; }}
  .copy-btn:hover {{ background:#f1f5f9; }}
  .footer {{ margin-top:28px; font-size:12px; color:#64748b; border-top:1px solid #e2e8f0; padding-top:12px; }}
</style>
</head>
<body>
<div class="container">
  <h1>CRISPR-Preflight Biological Dossier: {esc(report.gene_symbol)} <span style="font-weight:400;color:#64748b;">({esc(report.ensembl_id)})</span></h1>
  <div class="sub"><em>{esc(report.species_scientific_name)}</em> ({esc(report.species_common_name)}) &bull; {esc(report.gene_description)} &bull; {esc(report.assembly)} <code>{esc(report.chrom)}</code> ({'+' if report.strand >= 0 else '-'} strand) &bull; Primary transcript <code>{esc(tx.transcript_id)}</code> ({tx.num_exons} exons, {tx.protein_length_aa} aa)</div>

  <div class="verdict">
    <div style="font-size:12px;font-weight:700;text-transform:uppercase;">Recommended Perturbation Modality ({esc(mod.risk_level)} Trap Risk)</div>
    <div style="font-size:20px;font-weight:800;margin:4px 0;">{esc(mod.primary_modality)}</div>
    <div style="font-size:14px;">{esc(mod.summary_verdict)}</div>
  </div>

  <div class="kpi-grid">
    <div class="kpi"><div class="lbl">DepMap Cell Fitness (Chronos)</div><div class="val">{esc(chronos_metric)}</div>{ccle_sub}</div>
    <div class="kpi"><div class="lbl">gnomAD LoF Constraint (LOEUF)</div><div class="val">{esc(loeuf_metric)}</div></div>
    <div class="kpi"><div class="lbl">GTEx Baseline Expression</div><div class="val">{esc(tpm_metric)}</div></div>
    <div class="kpi"><div class="lbl">Coding Isoforms Audited</div><div class="val">{1 + len(report.alternative_isoforms)} isoform(s)</div></div>
  </div>

  <h2>1. Exon, Functional Domain &amp; Epitope Architecture</h2>
  {svg_diagram}
  {domain_key_html}

  <h2>2. Biological Trap Alerts</h2>
  {warnings_html}

  <h2>3. Exon-by-Exon Targetability &amp; Transparent Score Math</h2>
  <table>
    <thead>
      <tr><th>Exon</th><th>Score</th><th>Tier</th><th>AA Span</th><th>Length</th><th>Isoforms</th><th>NGG Sites</th><th>NMD Zone</th><th>Domain Overlap</th><th>CRISPick Coordinate</th></tr>
    </thead>
    <tbody>
      {''.join(exon_rows)}
    </tbody>
  </table>

  {guides_section}

  {crispri_section}

  {paralog_section}

  <h2>6. Wet-Lab Validation &amp; Orthogonal Chemical Controls</h2>
  <p><strong>Protein Assay:</strong> {esc(val.recommended_protein_assay)}</p>
  <p><strong>Western Antibody Selection:</strong> {esc(val.epitope_guidance_summary)}</p>
  <p><strong>RT-PCR Splicing Check:</strong> {esc(val.rt_qpcr_guidance)}</p>
  <h3 style="font-size:15px;margin-top:14px;">Orthogonal Chemical Probes &amp; Clinical Modulators</h3>
  {chem_items}

  <div class="footer">
    <div>Generated by <strong>crispr-preflight</strong> (developed by Mulham Najajreh) on {esc(report.generated_at)} &bull; Data sources: {esc(', '.join(report.data_sources_used))}. Please cite and respect each public resource's license terms.</div>
    <div style="margin-top:4px;">This work was supported in part by complimentary access to Google Antigravity and Gemini 4 Argon to accelerate scientific workflows.</div>
  </div>
</div>
<script id="preflight-data" type="application/json">{json_payload}</script>
<script>
function flashBtn(btn) {{
  const old = btn.dataset.origLabel || btn.textContent;
  btn.dataset.origLabel = old;
  btn.textContent = "Copied!";
  setTimeout(() => {{ btn.textContent = old; }}, 1400);
}}
function fallbackCopy(text, btn) {{
  const ta = document.createElement("textarea");
  ta.value = text;
  ta.style.position = "fixed";
  ta.style.opacity = "0";
  document.body.appendChild(ta);
  ta.select();
  try {{ document.execCommand("copy"); }} catch (e) {{}}
  document.body.removeChild(ta);
  flashBtn(btn);
}}
function copyCoord(btn) {{
  const text = btn.getAttribute("data-copy") || "";
  if (navigator.clipboard && navigator.clipboard.writeText) {{
    navigator.clipboard.writeText(text).then(() => flashBtn(btn)).catch(() => fallbackCopy(text, btn));
  }} else {{
    fallbackCopy(text, btn);
  }}
}}
</script>
</body>
</html>
"""
