"""Unit tests for modality decision rules (pan-essentiality, paralog buffering, silent baseline, validation)."""

from __future__ import annotations

import unittest

from crispr_preflight.engine.modality import classify_paralog_risks, evaluate_modality
from crispr_preflight.engine.validation import build_validation_plan
from crispr_preflight.models import (
    CellLineAlterationProfile,
    EssentialityProfile,
    ExonAudit,
    ExpressionProfile,
    IsoformSummary,
)


def _sample_isoform(num_exons: int = 8) -> IsoformSummary:
    return IsoformSummary(
        transcript_id="ENST000001",
        display_name="GENE-201",
        is_mane_select=True,
        is_mane_plus_clinical=False,
        is_canonical=True,
        tsl="1",
        cds_length_nt=1200,
        protein_length_aa=400,
        num_exons=num_exons,
        protein_id="ENSP000001",
    )


def _sample_exon(num: int = 2, tier: str = "RECOMMENDED", sym: bool = False) -> ExonAudit:
    return ExonAudit(
        exon_number=num,
        exon_id=f"ENSE{num}",
        chrom="chr1",
        genomic_start=1000,
        genomic_end=1140,
        strand=1,
        exon_length_bp=141,
        is_coding=True,
        coding_genomic_start=1000,
        coding_genomic_end=1140,
        coding_length_bp=141,
        aa_start=50,
        aa_end=97,
        is_symmetric=sym,
        nmd_status="COMPETENT",
        priority_score=95 if tier == "RECOMMENDED" else 55,
        tier=tier,
    )


class TestModalityAndValidationEngine(unittest.TestCase):
    def test_pan_essential_recommends_crispri_or_degron(self) -> None:
        ess = EssentialityProfile(
            gnomad_loeuf=0.18,
            depmap_total_lines=1100,
            depmap_essential_lines=1050,
            depmap_pan_essential_fraction=0.955,
            depmap_median_chronos=-1.82,
            matched_cell_line="A549",
            matched_cell_line_chronos=-1.91,
        )
        exp = ExpressionProfile(
            requested_tissue="lung",
            resolved_gtex_tissue="Lung",
            requested_tissue_tpm=45.2,
        )
        mod = evaluate_modality(
            gene_symbol="PLK1",
            primary_isoform=_sample_isoform(),
            alt_isoforms=[],
            essentiality=ess,
            expression=exp,
            paralogs=[],
            exons=[_sample_exon()],
        )
        self.assertEqual(mod.risk_level, "HIGH")
        self.assertIn("CRISPRi", mod.primary_modality)
        self.assertTrue(any(w.trap_id == "ESSENTIAL_CLONAL_ESCAPER" for w in mod.warnings))

    def test_coexpressed_paralog_recommends_dual_ko_control(self) -> None:
        raw_p = [
            {
                "gene_symbol": "SMARCA2",
                "ensembl_id": "ENSG00000080503",
                "percent_identity": 74.8,
                "target_percent_identity": 78.2,
                "ancestor_taxonomy": "Vertebrata",
            }
        ]
        risks = classify_paralog_risks(raw_p, {"SMARCA2": 18.4})
        self.assertEqual(risks[0].buffering_risk, "HIGH")

        mod = evaluate_modality(
            gene_symbol="SMARCA4",
            primary_isoform=_sample_isoform(),
            alt_isoforms=[],
            essentiality=EssentialityProfile(matched_cell_line="A549", matched_cell_line_chronos=-0.22),
            expression=ExpressionProfile("lung", "Lung", 32.0),
            paralogs=risks,
            exons=[_sample_exon()],
        )
        self.assertEqual(mod.risk_level, "HIGH")
        self.assertIn("SMARCA2", mod.primary_modality)
        self.assertTrue(any(w.trap_id == "PARALOG_BUFFERING" for w in mod.warnings))

    def test_silent_baseline_tissue_flags_expression_warning(self) -> None:
        mod = evaluate_modality(
            gene_symbol="ALB",
            primary_isoform=_sample_isoform(),
            alt_isoforms=[],
            essentiality=EssentialityProfile(),
            expression=ExpressionProfile(
                requested_tissue="brain",
                resolved_gtex_tissue="Brain_Cortex",
                requested_tissue_tpm=0.04,
                top_tissues=[("Liver", 48200.0)],
            ),
            paralogs=[],
            exons=[_sample_exon()],
        )
        self.assertEqual(mod.risk_level, "HIGH")
        self.assertTrue(any(w.trap_id == "SILENT_BASELINE_EXPRESSION" for w in mod.warnings))

    def test_antibody_epitope_avoids_reinitiation_and_symmetric_exons(self) -> None:
        e1 = ExonAudit(
            exon_number=1,
            exon_id="E1",
            chrom="chr1",
            genomic_start=100,
            genomic_end=250,
            strand=1,
            exon_length_bp=151,
            is_coding=True,
            aa_start=1,
            aa_end=50,
            reinitiation_met_codon=38,
            tier="CAUTION",
        )
        e2_sym = _sample_exon(num=2, tier="CAUTION", sym=True)
        e3_safe = ExonAudit(
            exon_number=3,
            exon_id="E3",
            chrom="chr1",
            genomic_start=2000,
            genomic_end=2900,
            strand=1,
            exon_length_bp=901,
            is_coding=True,
            aa_start=98,
            aa_end=400,
            is_symmetric=False,
            is_constitutive=True,
            tier="RECOMMENDED",
        )
        plan = build_validation_plan(
            gene_symbol="TEST",
            primary_isoform=_sample_isoform(3),
            exons=[e1, e2_sym, e3_safe],
            expression=ExpressionProfile(hpa_subcellular_locations=["Plasma membrane"]),
            chemical_tools=[],
        )
        self.assertIn("Flow Cytometry", plan.recommended_protein_assay)
        avoid_ranges = [(w.aa_start, w.aa_end) for w in plan.epitope_windows if w.status == "AVOID"]
        self.assertIn((1, 38), avoid_ranges)
        self.assertIn((50, 97), avoid_ranges)
        safe_ranges = [(w.aa_start, w.aa_end) for w in plan.epitope_windows if w.status == "SAFE"]
        self.assertIn((98, 400), safe_ranges)

    def test_ccle_existing_lof_and_high_amplification_traps(self) -> None:
        """Verify cBioPortal CCLE endogenous LoF (e.g. SMARCA4 in A549) and GISTIC +2 amplification traps."""
        lof_mod = evaluate_modality(
            gene_symbol="SMARCA4",
            primary_isoform=_sample_isoform(),
            alt_isoforms=[],
            essentiality=EssentialityProfile(matched_cell_line="A549", matched_cell_line_chronos=-0.12),
            expression=ExpressionProfile("lung", "Lung", 25.0),
            paralogs=[],
            exons=[_sample_exon()],
            cell_line_genetics=CellLineAlterationProfile(
                requested_cell_line="A549",
                matched_ccle_sample="A549_LUNG",
                mutations=["p.Q729Cfs*4 (Frame_Shift_Del, VAF 100%)"],
                has_truncating_lof=True,
                cna_gistic=-1,
                cna_label="Hemizygous Loss (GISTIC -1)",
            ),
        )
        self.assertEqual(lof_mod.risk_level, "HIGH")
        self.assertTrue(any(w.trap_id == "CELL_LINE_EXISTING_LOF" for w in lof_mod.warnings))

        amp_mod = evaluate_modality(
            gene_symbol="MYC",
            primary_isoform=_sample_isoform(),
            alt_isoforms=[],
            essentiality=EssentialityProfile(matched_cell_line="MCF7", matched_cell_line_chronos=-0.65),
            expression=ExpressionProfile("breast", "Breast_Mammary_Tissue", 80.0),
            paralogs=[],
            exons=[_sample_exon()],
            cell_line_genetics=CellLineAlterationProfile(
                requested_cell_line="MCF7",
                matched_ccle_sample="MCF7_BREAST",
                mutations=[],
                has_truncating_lof=False,
                cna_gistic=2,
                cna_label="High-Level Amplification (GISTIC +2)",
            ),
        )
        self.assertEqual(amp_mod.risk_level, "HIGH")
        self.assertIn("CRISPRi", amp_mod.primary_modality)
        self.assertTrue(any(w.trap_id == "CELL_LINE_HIGH_AMPLIFICATION" for w in amp_mod.warnings))


if __name__ == "__main__":
    unittest.main()
