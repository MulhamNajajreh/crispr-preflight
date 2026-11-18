"""Perturbation modality decision engine (Cas9 KO vs. CRISPRi/dTAG vs. Dual-Paralog KO).

Synthesizes DepMap Chronos cell fitness, gnomAD LOEUF constraint, GTEx/HPA baseline
expression, Ensembl Compara paralog buffering, and intronless/isoform architecture.
"""

from __future__ import annotations

from typing import Optional

from crispr_preflight.models import (
    CellLineAlterationProfile,
    EssentialityProfile,
    ExonAudit,
    ExpressionProfile,
    IsoformSummary,
    ModalityRecommendation,
    ParalogRisk,
    TrapWarning,
)


def classify_paralog_risks(
    raw_paralogs: list[dict],
    paralog_tissue_tpms: dict[str, Optional[float]],
) -> list[ParalogRisk]:
    """Assign HIGH / MODERATE / LOW buffering risk based on sequence identity and tissue co-expression."""
    risks: list[ParalogRisk] = []
    for p in raw_paralogs:
        sym = str(p["gene_symbol"])
        q_id = float(p["percent_identity"])
        t_id = float(p["target_percent_identity"])
        max_id = max(q_id, t_id)
        tpm = paralog_tissue_tpms.get(sym)

        # De Kegel et al. (Cell Systems 2021): high buffering above ~55-60% identity when co-expressed
        if max_id >= 55.0 and (tpm is None or tpm >= 2.0):
            level = "HIGH" if (tpm is not None and tpm >= 2.0) or max_id >= 70.0 else "MODERATE"
        elif max_id >= 40.0 and (tpm is None or tpm >= 1.5):
            level = "MODERATE"
        else:
            level = "LOW"

        risks.append(
            ParalogRisk(
                gene_symbol=sym,
                ensembl_id=str(p["ensembl_id"]),
                percent_identity=q_id,
                target_percent_identity=t_id,
                ancestor_taxonomy=p.get("ancestor_taxonomy"),
                target_tissue_tpm=tpm,
                buffering_risk=level,
            )
        )
    risk_order = {"HIGH": 0, "MODERATE": 1, "LOW": 2}
    risks.sort(key=lambda r: (risk_order.get(r.buffering_risk, 9), -max(r.percent_identity, r.target_percent_identity)))
    return risks


def evaluate_modality(
    gene_symbol: str,
    primary_isoform: IsoformSummary,
    alt_isoforms: list[IsoformSummary],
    essentiality: EssentialityProfile,
    expression: ExpressionProfile,
    paralogs: list[ParalogRisk],
    exons: list[ExonAudit],
    cell_line_genetics: Optional[CellLineAlterationProfile] = None,
    species_id: str = "human",
    human_ortholog_symbol: Optional[str] = None,
    nuclease_id: str = "spcas9",
) -> ModalityRecommendation:
    """Produce a transparent, rule-based experimental modality verdict and trap warning list."""
    warnings: list[TrapWarning] = []
    coding_exons = [e for e in exons if e.is_coding]
    recommended_exons = [e for e in coding_exons if e.tier == "RECOMMENDED"]

    # --- 0. Endogenous Cell-Line Mutation & Copy-Number Trap (CCLE / DepMap) ---
    cell_existing_lof = False
    cell_high_amp = False
    cell_label = ""
    if cell_line_genetics and cell_line_genetics.matched_ccle_sample:
        cell_label = cell_line_genetics.matched_ccle_sample.split("_", 1)[0]
        if cell_line_genetics.has_truncating_lof or cell_line_genetics.cna_gistic == -2:
            cell_existing_lof = True
            alt_bits = list(cell_line_genetics.mutations)
            if cell_line_genetics.cna_gistic == -2:
                alt_bits.append("Homozygous Deletion (GISTIC -2)")
            alt_str = "; ".join(alt_bits) or "Loss-of-function alteration"
            warnings.append(
                TrapWarning(
                    trap_id="CELL_LINE_EXISTING_LOF",
                    title=f"Already LoF-Mutated or Deleted in {cell_label}: {alt_str}",
                    severity="HIGH",
                    explanation=(
                        f"In CCLE/DepMap ({cell_line_genetics.matched_ccle_sample}), {gene_symbol} already carries "
                        f"an endogenous loss-of-function alteration ({alt_str}). Attempting a 'knockout' in a line "
                        "that lacks wild-type functional protein will yield false-negative or confounded results."
                    ),
                    bench_action=(
                        f"Confirm full-length wild-type {gene_symbol} protein by Western blot in {cell_label} before "
                        f"perturbation, choose a {gene_symbol}-wild-type cell line, or perform cDNA wild-type add-back."
                    ),
                )
            )
        elif cell_line_genetics.mutations:
            mut_str = "; ".join(cell_line_genetics.mutations[:3])
            warnings.append(
                TrapWarning(
                    trap_id="CELL_LINE_SOMATIC_MUTATION",
                    title=f"Endogenous Coding Variant in {cell_label}: {mut_str}",
                    severity="MODERATE",
                    explanation=(
                        f"{cell_label} ({cell_line_genetics.matched_ccle_sample}) harbors endogenous {gene_symbol} "
                        f"coding mutation(s): {mut_str}. If an sgRNA overlaps this mutated codon, wild-type reference "
                        "spacers will mismatch the mutant allele; oncogenic gain-of-function drivers also create oncogene addiction."
                    ),
                    bench_action=(
                        "Avoid placing sgRNAs directly over the mutated residue unless intentionally designing an allele-specific guide."
                    ),
                )
            )

        if cell_line_genetics.cna_gistic == 2:
            cell_high_amp = True
            warnings.append(
                TrapWarning(
                    trap_id="CELL_LINE_HIGH_AMPLIFICATION",
                    title=f"High-Level Focal Amplification in {cell_label} (GISTIC +2): Cas9 Multi-DSB Toxicity Trap",
                    severity="HIGH",
                    explanation=(
                        f"{gene_symbol} sits in a high-copy-number amplicon in {cell_label}. Cutting an amplified locus "
                        "with Cas9 nuclease generates multiple simultaneous DNA double-strand breaks, triggering gene-independent "
                        "p53/G2-M cell-cycle arrest and incomplete multi-allele editing (Aguirre et al., 2016; Meyers et al., 2017)."
                    ),
                    bench_action=(
                        "Use dCas9-KRAB (CRISPRi) promoter repression instead of Cas9 nuclease to avoid copy-number DSB toxicity."
                    ),
                )
            )

    # --- 1. Baseline Expression Check (Silent Target Trap) ---
    low_expression = False
    if expression.requested_tissue and expression.requested_tissue_tpm is not None:
        if expression.requested_tissue_tpm < 1.0:
            low_expression = True
            top_str = (
                ", ".join(f"{t} ({v} TPM)" for t, v in expression.top_tissues[:3])
                if expression.top_tissues
                else "none >= 1 TPM in GTEx"
            )
            warnings.append(
                TrapWarning(
                    trap_id="SILENT_BASELINE_EXPRESSION",
                    title=f"Very Low Baseline Expression in '{expression.resolved_gtex_tissue}' ({expression.requested_tissue_tpm} TPM)",
                    severity="HIGH",
                    explanation=(
                        f"{gene_symbol} has <1.0 median TPM in normal '{expression.resolved_gtex_tissue}'. "
                        f"Knocking out an unexpressed gene yields a silent null unless expression is stimulus-induced "
                        f"or up-regulated in your specific tumor line. (Top GTEx tissues: {top_str})."
                    ),
                    bench_action=(
                        "Verify baseline expression in your exact cell line by RT-qPCR/RNA-seq before ordering KO guides, "
                        "or use CRISPRa if studying gain-of-function."
                    ),
                )
            )
    elif expression.top_tissues and expression.top_tissues[0][1] < 1.0:
        low_expression = True
        warnings.append(
            TrapWarning(
                trap_id="RESTRICTED_EXPRESSION",
                title="Restricted / Near-Zero Adult Tissue Expression (<1.0 TPM across GTEx)",
                severity="MODERATE",
                explanation=(
                    f"{gene_symbol} peaks below 1.0 TPM across all 54 adult GTEx tissues "
                    "(may be developmental, inducible, or lineage-restricted)."
                ),
                bench_action="Confirm mRNA/protein presence in your target cells prior to perturbation.",
            )
        )

    # --- 2. Essentiality & Clonal In-Frame Escaper Trap ---
    cell_chronos = essentiality.matched_cell_line_chronos
    lin_chronos = essentiality.matched_lineage_median_chronos
    pan_frac = essentiality.depmap_pan_essential_fraction
    loeuf = essentiality.gnomad_loeuf

    severe_essential = False
    moderate_essential = False

    if cell_chronos is not None and cell_chronos <= -0.85:
        severe_essential = True
        warnings.append(
            TrapWarning(
                trap_id="ESSENTIAL_CLONAL_ESCAPER",
                title=(
                    f"Strongly Essential in {essentiality.matched_cell_line} "
                    f"(DepMap Chronos = {cell_chronos:+.2f})"
                ),
                severity="HIGH",
                explanation=(
                    f"Cas9 induces ~1/3 in-frame indels (-3 bp, -6 bp). Because true biallelic loss of {gene_symbol} "
                    f"is lethal/anti-proliferative in {essentiality.matched_cell_line}, single-cell clonal KO expansion "
                    "strongly selects for surviving in-frame hypomorphic escapers."
                ),
                bench_action=(
                    "Avoid 3-week clonal Cas9 KO isolation. Use CRISPRi knockdown, an inducible degron (dTAG/AID), "
                    "or assay pooled Cas9 populations early (Day 4–7 post-transduction) while targeting a catalytic domain."
                ),
            )
        )
    elif (pan_frac is not None and pan_frac >= 0.75) or bool(essentiality.is_core_essential_ot):
        severe_essential = True
        pct = int((pan_frac or 0.8) * 100)
        orth_tag = (
            f" [via 1:1 human ortholog {human_ortholog_symbol}]"
            if (species_id != "human" and human_ortholog_symbol)
            else ""
        )
        warnings.append(
            TrapWarning(
                trap_id="PAN_ESSENTIAL_ESCAPER",
                title=f"Core Pan-Essential Fitness Gene{orth_tag} ({pct}% DepMap lines <= -1.0)",
                severity="HIGH",
                explanation=(
                    f"{gene_symbol}{orth_tag} behaves as a deeply conserved core cell-fitness gene "
                    f"(median Chronos = {essentiality.depmap_median_chronos}). Prolonged culture selection or single-clone "
                    "isolation strongly enriches for in-frame (-3 bp) hypomorphic escapers or unedited cells; "
                    "homozygous whole-organism KO is frequently embryonic/larval lethal."
                ),
                bench_action=(
                    "Use CRISPRi/conditional Cre-lox/RNAi knockdown, acute degron tagging (dTAG/AID), or assay F0/pooled "
                    "populations early while targeting a conserved catalytic domain."
                ),
            )
        )
    elif (
        (cell_chronos is not None and -0.85 < cell_chronos <= -0.5)
        or (lin_chronos is not None and lin_chronos <= -0.65)
        or (loeuf is not None and loeuf < 0.35)
    ):
        moderate_essential = True
        detail_parts: list[str] = []
        if cell_chronos is not None:
            detail_parts.append(f"{essentiality.matched_cell_line} Chronos = {cell_chronos:+.2f}")
        elif lin_chronos is not None:
            detail_parts.append(f"{essentiality.matched_lineage} median Chronos = {lin_chronos:+.2f}")
        if loeuf is not None:
            detail_parts.append(f"gnomAD LOEUF = {loeuf:.2f}")
        warnings.append(
            TrapWarning(
                trap_id="CONSTRAINT_FITNESS_RISK",
                title=f"Evolutionary Constraint / Moderate Fitness Cost ({', '.join(detail_parts)})",
                severity="MODERATE",
                explanation=(
                    f"{gene_symbol} is intolerant to loss-of-function in humans or confers a moderate proliferation "
                    "penalty in culture, meaning wild-type or in-frame indel alleles can outgrow frameshift clones over passages."
                ),
                bench_action=(
                    "Target a conserved functional domain (so -3 bp indels are also deleterious) and quantify allele "
                    "editing % by ICE/amplicon-seq at both early (Day 4) and assay timepoints."
                ),
            )
        )

    # --- 3. Paralog Buffering / Functional Redundancy Trap ---
    high_paralogs = [p for p in paralogs if p.buffering_risk == "HIGH"]
    mod_paralogs = [p for p in paralogs if p.buffering_risk == "MODERATE"]

    if high_paralogs:
        p_desc = ", ".join(
            f"{p.gene_symbol} ({p.percent_identity:.0f}% id"
            + (f", {p.target_tissue_tpm} TPM)" if p.target_tissue_tpm is not None else ")")
            for p in high_paralogs[:3]
        )
        warnings.append(
            TrapWarning(
                trap_id="PARALOG_BUFFERING",
                title=f"High Paralog Compensation Risk: {p_desc}",
                severity="HIGH",
                explanation=(
                    f"Close human paralog(s) share high sequence identity and are expressed in the target context. "
                    f"Single-gene {gene_symbol} knockout frequently yields muted or false-negative phenotypes due to paralog buffering (De Kegel et al., 2021)."
                ),
                bench_action=(
                    f"Include a dual-guide digenic KO control co-targeting {gene_symbol} + {high_paralogs[0].gene_symbol} "
                    "alongside single-gene perturbations."
                ),
            )
        )
    elif mod_paralogs:
        p_desc = ", ".join(f"{p.gene_symbol} ({p.percent_identity:.0f}% id)" for p in mod_paralogs[:3])
        warnings.append(
            TrapWarning(
                trap_id="PARALOG_MODERATE",
                title=f"Related Human Paralog(s) Present: {p_desc}",
                severity="MODERATE",
                explanation=(
                    f"Moderate-identity paralog(s) exist for {gene_symbol} and can partially buffer loss-of-function depending on cell-line expression."
                ),
                bench_action=(
                    f"Check RNA expression of {mod_paralogs[0].gene_symbol} in your model line; consider dual perturbation if single KO shows no phenotype."
                ),
            )
        )

    # --- 4. Intronless / Single-Exon Global NMD Escape Trap ---
    is_single_exon = primary_isoform.num_exons == 1
    if is_single_exon:
        warnings.append(
            TrapWarning(
                trap_id="SINGLE_EXON_NMD_ESCAPE",
                title=f"Single-Exon Transcript ({primary_isoform.transcript_id}): Global EJC-NMD Escape",
                severity="HIGH",
                explanation=(
                    f"{gene_symbol} has only 1 exon (0 splice junctions -> 0 Exon Junction Complexes). "
                    "Frameshift indels will NEVER trigger EJC-dependent Nonsense-Mediated mRNA Decay; "
                    "mutated mRNA remains stable and translates a truncated polypeptide."
                ),
                bench_action=(
                    "Do NOT cut near the middle/3' CDS expecting mRNA decay. Either (1) use CRISPRi to silence transcription, "
                    "(2) use paired dual-guides to excise the entire ORF, or (3) cut near the 5' start inside a critical folded domain."
                ),
            )
        )

    # --- 5. Alternative Isoform / Symmetric Exon Complexity ---
    skipped_exons = [e for e in coding_exons if not e.is_constitutive]
    symmetric_exons = [e for e in coding_exons if e.is_symmetric]
    if alt_isoforms and skipped_exons:
        warnings.append(
            TrapWarning(
                trap_id="ISOFORM_SPLICING_TRAP",
                title=f"{len(skipped_exons)}/{len(coding_exons)} Coding Exon(s) Spliced Out in Alternative Isoforms",
                severity="MODERATE",
                explanation=(
                    f"{gene_symbol} expresses {1 + len(alt_isoforms)} well-supported coding isoforms. "
                    f"Exon(s) {', '.join(str(e.exon_number) for e in skipped_exons[:5])} are absent in at least one coding isoform."
                ),
                bench_action="Restrict sgRNA design strictly to constitutive exons shared across all expressed coding isoforms.",
            )
        )

    if symmetric_exons:
        warnings.append(
            TrapWarning(
                trap_id="SYMMETRIC_EXON_SKIPPING",
                title=f"{len(symmetric_exons)} Symmetric Internal Exon(s) (Length Divisible by 3)",
                severity="INFO" if recommended_exons else "MODERATE",
                explanation=(
                    f"Exon(s) {', '.join(str(e.exon_number) for e in symmetric_exons[:6])} have coding lengths divisible by 3. "
                    "Indels disrupting Exonic Splicing Enhancers (ESEs) can trigger clean in-frame exon skipping (Tuladhar et al., 2019; Smits et al., 2019)."
                ),
                bench_action=(
                    "Prefer asymmetric exons or exons encoding catalytic residues, and place RT-PCR validation primers spanning flanking exons to catch exon skipping."
                ),
            )
        )

    if species_id == "zebrafish":
        warnings.append(
            TrapWarning(
                trap_id="SPECIES_ZEBRAFISH_NOTE",
                title="Zebrafish (GRCz11) Bench Note: Teleost Ohnologs & F0 Crispant Design",
                severity="INFO",
                explanation=(
                    "Teleost fishes underwent an ancestral whole-genome duplication (~25% of genes retain 'a'/'b' ohnolog pairs) "
                    "and exhibit high polymorphism across AB/TU/TL lab strains."
                ),
                bench_action=(
                    "For F0 crispants, co-inject 2–3 non-overlapping Cas9 RNPs across top asymmetric exons (and co-target any high-identity 'a'/'b' paralog above); "
                    "sequence your wild-type clutch with the TIDE/ICE PCR primers before microinjection."
                ),
            )
        )
    elif species_id == "fly":
        warnings.append(
            TrapWarning(
                trap_id="SPECIES_FLY_NOTE",
                title="Drosophila (BDGP6.46) Bench Note: U6 Vector Overhangs & Strain SNPs",
                severity="INFO",
                explanation=(
                    "For Drosophila germline pCFD3 (dU6:3) cloning, BbsI overhangs use 5'-GTCG-3' (Fwd) / 5'-AAAC-3' (Rev) "
                    "instead of mammalian 5'-CACC-3'; for S2R+ RNP or embryo injection, order the bare 20-nt crRNA/IVT template."
                ),
                bench_action=(
                    "Sequence the target locus in your exact injection stock (e.g. nos-Cas9 / vasa-Cas9) using the paired PCR primers to rule out balancer/strain SNPs."
                ),
            )
        )
    elif species_id == "worm":
        warnings.append(
            TrapWarning(
                trap_id="SPECIES_WORM_NOTE",
                title="C. elegans (WBcel235) Bench Note: Gonad RNP Delivery & dpy-10 Co-CRISPR",
                severity="INFO",
                explanation=(
                    "Germline syncytial gonad microinjection in C. elegans is typically performed with synthetic RNP complexes "
                    "rather than mammalian U6 lentiviral vectors."
                ),
                bench_action=(
                    "Pair your top crRNA spacer with a dpy-10 or pha-1 co-CRISPR marker and genotype F1 rollers/dumpies using the paired PCR primers."
                ),
            )
        )

    # --- 6. Nuclease-Specific Biochemical & Cloning Notes ---
    if nuclease_id == "spcas9_hf":
        warnings.append(
            TrapWarning(
                trap_id="NUCLEASE_HIFI_5G_RULE",
                title="High-Fidelity Cas9 (eSpCas9 / SpCas9-HF1 / HiFi): Native 5'-G Rule",
                severity="INFO",
                explanation=(
                    "Engineered eSpCas9(1.1) and SpCas9-HF1 lose 50–90% activity under standard U6 promoters when an unmatched "
                    "5'-G is appended to a 20-nt spacer (21-nt gN20; Kim et al., Nat Commun 2017; Kulcsar et al., 2017)."
                ),
                bench_action=(
                    "Pick native 5'-GN19 guides ranked first below for U6 plasmid cloning, or deliver non-G spacers as synthetic "
                    "IDT Alt-R HiFi Cas9 RNPs (or via a tRNA/Hammerhead self-cleaving sgRNA cassette)."
                ),
            )
        )
    elif nuclease_id == "spg":
        warnings.append(
            TrapWarning(
                trap_id="NUCLEASE_SPG_OFFTARGET",
                title="PAM-Relaxed SpG Cas9 (3'-NGN): Expanded Off-Target Check Required",
                severity="INFO",
                explanation=(
                    "SpG recognizes 3'-NGN PAMs (4x denser than canonical NGG; Walton et al., Science 2020), rescuing compact "
                    "catalytic exons while increasing candidate genome-wide off-target loci."
                ),
                bench_action=(
                    "Run Cas-OFFinder with PAM set to 'NGN' (0–2 mismatches) before ordering oligos, or pair with HiFi-SpG / RNP delivery."
                ),
            )
        )
    elif nuclease_id == "sacas9":
        warnings.append(
            TrapWarning(
                trap_id="NUCLEASE_SACAS9_AAV",
                title="Compact SaCas9 (21-nt + 3'-NNGRRT): Single-AAV Packaging & BsaI Cloning",
                severity="INFO",
                explanation=(
                    "S. aureus Cas9 (~3.16 kb; Ran et al., Nature 2015) fits with its U6-sgRNA cassette inside a single all-in-one "
                    "AAV vector (~4.7 kb limit, e.g. pX601 Addgene #61591) using 21-nt spacers and BsaI Golden Gate overhangs."
                ),
                bench_action=(
                    "Digest pX601 with BsaI (not BsmBI/BbsI) when cloning the CACC/AAAC oligos below, and verify off-targets in Cas-OFFinder with PAM = 'NNGRRT'."
                ),
            )
        )
    elif nuclease_id == "cas12a":
        ectotherm = species_id in ("zebrafish", "fly", "worm")
        warnings.append(
            TrapWarning(
                trap_id="NUCLEASE_CAS12A_NOTE",
                title=(
                    "Cas12a / Cpf1 (5'-TTTV + 23-nt): Use LbCas12a-Ultra / enAsCas12a Below 30 °C"
                    if ectotherm
                    else "Cas12a / Cpf1 (5'-TTTV + 23-nt): Staggered 5-nt Overhangs & crRNA Multiplexing"
                ),
                severity="MODERATE" if ectotherm else "INFO",
                explanation=(
                    "Wild-type AsCas12a has sharply reduced cleavage activity below 28–30 °C in ectothermic embryos/animals "
                    "(Moreno-Mateos et al., Nat Commun 2017), whereas LbCas12a (LbCas12a-Ultra) and enAsCas12a remain highly active."
                    if ectotherm
                    else "Cas12a recognizes a 5' T-rich PAM (5'-TTTV) upstream of the 23-nt protospacer, creates staggered 5-nt overhangs "
                    "that favor frame-disrupting microdeletions, and processes multiplexed crRNA arrays from a single transcript."
                ),
                bench_action=(
                    "Use LbCas12a-Ultra or enAsCas12a RNP (with 34 °C brief heat-shock in zebrafish if applicable); cloning oligos below use BsmBI crRNA direct-repeat overhangs (5'-AGAT / 5'-AAAA)."
                    if ectotherm
                    else "For U6 crRNA plasmids (pY010 / pY016 / pRDA_052), anneal the BsmBI 5'-AGAT / 5'-AAAA oligos below; for synthetic IDT Alt-R crRNA RNP, use the 21–23 nt protospacer."
                ),
            )
        )

    # --- Synthesize Top-Level Modality Verdict ---
    if cell_existing_lof:
        primary_mod = f"Switch Cell Line or Wild-Type Rescue ({gene_symbol} Already LoF in {cell_label})"
        risk_level = "HIGH"
        verdict = (
            f"{gene_symbol} already carries a truncating mutation or homozygous deletion in {cell_label} (CCLE). "
            "Verify wild-type full-length protein expression before attempting knockout, or switch to a wild-type model line."
        )
    elif cell_high_amp:
        primary_mod = f"CRISPRi Knockdown ({gene_symbol} Amplified in {cell_label} — Avoid Cas9 DSB Toxicity)"
        risk_level = "HIGH"
        verdict = (
            f"{gene_symbol} is high-level amplified (GISTIC +2) in {cell_label}; Cas9 nuclease cutting causes "
            "multi-break DNA damage toxicity. Use dCas9-KRAB (CRISPRi) promoter guides instead."
        )
    elif low_expression and expression.requested_tissue:
        primary_mod = "Verify Cell-Line Expression First (or CRISPRa)"
        risk_level = "HIGH"
        verdict = (
            f"Baseline {gene_symbol} expression in '{expression.resolved_gtex_tissue}' is <1 TPM. "
            "Confirm expression in your specific cell line/condition before running a loss-of-function KO."
        )
    elif severe_essential:
        primary_mod = "CRISPRi / Inducible Degron (dTAG) / Acute Pooled Cas9"
        risk_level = "HIGH"
        verdict = (
            f"{gene_symbol} is strongly fitness-essential; clonal Cas9 KO selection risks enriching for -3 bp in-frame escapers. "
            "Prefer tunable knockdown (CRISPRi), acute degradation (dTAG), or early pooled domain-targeted Cas9 assays."
        )
    elif is_single_exon:
        primary_mod = "CRISPRi or Dual-Guide Whole-ORF Excision (or 5' Domain KO)"
        risk_level = "HIGH"
        verdict = (
            f"{gene_symbol} is a single-exon transcript where frameshifts escape NMD. "
            "Use CRISPRi, dual-guide deletion, or target an essential N-terminal domain."
        )
    elif high_paralogs:
        top_par = high_paralogs[0].gene_symbol
        primary_mod = f"Single Cas9 KO + Dual-Guide ({gene_symbol} + {top_par}) Control"
        risk_level = "HIGH"
        verdict = (
            f"{gene_symbol} is safe for Cas9 KO ({len(recommended_exons)} recommended exon(s)), "
            f"but co-expressed paralog {top_par} poses high functional buffering risk—include a dual-KO arm."
        )
    elif moderate_essential or mod_paralogs or not recommended_exons:
        primary_mod = "Domain-Targeted Cas9 KO (with Allele & Isoform QC)"
        risk_level = "MODERATE"
        verdict = (
            f"Cas9 KO is feasible ({len(recommended_exons)} top-tier exon(s)), "
            "with moderate constraint, isoform, or paralog nuances flagged below."
        )
    else:
        primary_mod = "Standard Domain-Targeted Cas9 KO"
        risk_level = "LOW"
        verdict = (
            f"{gene_symbol} is a clean Cas9 KO candidate with {len(recommended_exons)} high-confidence "
            "constitutive, NMD-competent coding exon(s) and low escaper/paralog risk."
        )

    return ModalityRecommendation(
        primary_modality=primary_mod,
        risk_level=risk_level,
        summary_verdict=verdict,
        warnings=warnings,
    )
