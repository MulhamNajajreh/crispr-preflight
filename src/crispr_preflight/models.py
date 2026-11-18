"""Data models for crispr-preflight reports, exons, isoforms, and biological risk assessments.

Uses pure Python standard-library dataclasses so the package installs cleanly
on any Python >=3.9 environment with zero compiled C/Rust dependencies.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Optional


@dataclass
class IsoformSummary:
    """Metadata summary for a protein-coding Ensembl transcript isoform."""

    transcript_id: str
    display_name: str
    is_mane_select: bool
    is_mane_plus_clinical: bool
    is_canonical: bool
    tsl: Optional[str]
    cds_length_nt: int
    protein_length_aa: int
    num_exons: int
    protein_id: Optional[str] = None


@dataclass
class ProteinDomain:
    """Functional or structural domain annotation mapped onto protein coordinates."""

    source_db: str
    accession: str
    description: str
    aa_start: int
    aa_end: int


@dataclass
class GenotypingPrimerPair:
    """Genomic PCR + Sanger sequencing primer pair optimized for TIDE / Synthego ICE indel decomposition."""

    fwd_seq: str
    rev_seq: str
    fwd_tm: float
    rev_tm: float
    amplicon_bp: int
    cut_offset_from_fwd: int
    cut_offset_from_rev: int
    recommended_seq_primer: str  # "Forward" or "Reverse"
    genomic_region: str


@dataclass
class GuideCandidate:
    """Candidate 20-nt SpCas9 sgRNA spacer for coding KO or TSS-proximal CRISPRi.

    On-target sequence/context heuristics are evaluated locally; genome-wide
    off-target alignment (0-3 mismatches across GRCh38) should be verified in
    Broad CRISPick or Cas-OFFinder prior to ordering oligos.
    """

    guide_id: str
    exon_number: int
    spacer_20nt: str
    pam: str
    genomic_strand: str  # "+" or "-" on GRCh38
    transcript_orientation: str  # "Sense" or "Antisense" relative to mRNA
    chrom: str
    genomic_cut_pos: int
    aa_cut: int
    gc_percent: int
    overlaps_domain: Optional[str] = None
    sequence_score: int = 0
    quality_notes: list[str] = field(default_factory=list)
    cloning_fwd_oligo: str = ""
    cloning_rev_oligo: str = ""
    mode: str = "KO"  # "KO" or "CRISPRi"
    tss_offset_bp: Optional[int] = None
    genomic_start_23bp: Optional[int] = None
    genomic_end_23bp: Optional[int] = None
    common_snp_note: Optional[str] = None
    genotyping_primers: Optional[GenotypingPrimerPair] = None
    nuclease: str = "SpCas9"
    spacer_length: int = 20
    pam_orientation: str = "3prime"  # "3prime" ([spacer][PAM]) or "5prime" ([PAM][spacer])
    cloning_backbone_hint: str = "U6 BsmBI/BbsI (lentiCRISPRv2 / PX459)"


@dataclass
class ExonAudit:
    """Biological and geometric audit for a single exon in 5' -> 3' transcript order."""

    exon_number: int
    exon_id: str
    chrom: str
    genomic_start: int
    genomic_end: int
    strand: int
    exon_length_bp: int
    is_coding: bool
    coding_genomic_start: Optional[int] = None
    coding_genomic_end: Optional[int] = None
    coding_length_bp: int = 0
    cds_offset_start: Optional[int] = None
    cds_offset_end: Optional[int] = None
    aa_start: Optional[int] = None
    aa_end: Optional[int] = None
    is_constitutive: bool = True
    isoform_coverage_fraction: float = 1.0
    skipped_in_isoforms: list[str] = field(default_factory=list)
    is_symmetric: bool = False
    nmd_status: str = "NON_CODING_UTR"
    nmd_escape_fraction: float = 0.0
    reinitiation_met_codon: Optional[int] = None
    overlapping_domains: list[str] = field(default_factory=list)
    spcas9_pam_count: Optional[int] = None
    priority_score: int = 0
    tier: str = "UTR_NONCODING"
    score_breakdown: list[str] = field(default_factory=list)
    crispick_target_region: Optional[str] = None
    candidate_guides: list[GuideCandidate] = field(default_factory=list)


@dataclass
class ParalogRisk:
    """Assessment of functional buffering risk from a human paralog."""

    gene_symbol: str
    ensembl_id: str
    percent_identity: float
    target_percent_identity: float
    ancestor_taxonomy: Optional[str] = None
    target_tissue_tpm: Optional[float] = None
    buffering_risk: str = "LOW"


@dataclass
class EssentialityProfile:
    """Combined organism-level (gnomAD) and cell-level (DepMap Chronos) essentiality."""

    gnomad_pli: Optional[float] = None
    gnomad_loeuf: Optional[float] = None
    gnomad_oe_lof: Optional[float] = None
    is_core_essential_ot: Optional[bool] = None
    depmap_total_lines: int = 0
    depmap_essential_lines: int = 0
    depmap_pan_essential_fraction: Optional[float] = None
    depmap_median_chronos: Optional[float] = None
    matched_cell_line: Optional[str] = None
    matched_cell_line_lineage: Optional[str] = None
    matched_cell_line_chronos: Optional[float] = None
    matched_lineage: Optional[str] = None
    matched_lineage_median_chronos: Optional[float] = None
    matched_lineage_cell_count: int = 0


@dataclass
class ExpressionProfile:
    """Baseline mRNA (GTEx) and protein/localization (HPA) expression context."""

    requested_tissue: Optional[str] = None
    resolved_gtex_tissue: Optional[str] = None
    requested_tissue_tpm: Optional[float] = None
    top_tissues: list[tuple[str, float]] = field(default_factory=list)
    hpa_subcellular_locations: list[str] = field(default_factory=list)
    hpa_tissue_ihc_level: Optional[str] = None
    hpa_rna_specificity: Optional[str] = None


@dataclass
class ChemicalTool:
    """Orthogonal small-molecule chemical probe or clinical modulator."""

    name: str
    identifier: str
    tool_type: str
    mechanism_of_action: str
    quality_or_stage: str
    has_negative_control: bool = False
    source_url: Optional[str] = None


@dataclass
class TrapWarning:
    """Structured biological failure-mode alert with actionable bench guidance."""

    trap_id: str
    title: str
    severity: str  # "HIGH", "MODERATE", "INFO"
    explanation: str
    bench_action: str


@dataclass
class ModalityRecommendation:
    """Top-level perturbation strategy recommendation before designing sgRNAs."""

    primary_modality: str
    risk_level: str  # "HIGH", "MODERATE", "LOW"
    summary_verdict: str
    warnings: list[TrapWarning] = field(default_factory=list)


@dataclass
class EpitopeWindow:
    """Amino-acid interval annotated as safe or unsafe for Western blot antibody binding."""

    aa_start: int
    aa_end: int
    status: str  # "SAFE" or "AVOID"
    reason: str


@dataclass
class ValidationPlan:
    """Wet-lab readout, antibody epitope guidance, and orthogonal chemical controls."""

    recommended_protein_assay: str
    epitope_guidance_summary: str
    epitope_windows: list[EpitopeWindow] = field(default_factory=list)
    rt_qpcr_guidance: str = ""
    orthogonal_chemical_tools: list[ChemicalTool] = field(default_factory=list)


@dataclass
class CellLineAlterationProfile:
    """Endogenous somatic mutation and copy-number profile of the target gene in a CCLE cancer cell line."""

    requested_cell_line: Optional[str] = None
    matched_ccle_sample: Optional[str] = None
    mutations: list[str] = field(default_factory=list)
    has_truncating_lof: bool = False
    cna_gistic: Optional[int] = None  # -2 homdel, -1 hetloss, 0 neutral, 1 gain, 2 high-amp
    cna_label: Optional[str] = None


@dataclass
class PreflightReport:
    """Complete pre-experiment biological flight check dossier for a target gene."""

    gene_symbol: str
    ensembl_id: str
    gene_description: str
    chrom: str
    strand: int
    assembly: str
    primary_transcript: IsoformSummary
    alternative_isoforms: list[IsoformSummary]
    essentiality: EssentialityProfile
    expression: ExpressionProfile
    paralogs: list[ParalogRisk]
    modality: ModalityRecommendation
    domains: list[ProteinDomain]
    exons: list[ExonAudit]
    validation: ValidationPlan
    data_sources_used: list[str]
    generated_at: str
    recommended_guides: list[GuideCandidate] = field(default_factory=list)
    crispri_guides: list[GuideCandidate] = field(default_factory=list)
    tss_genomic_pos: Optional[int] = None
    cell_line_genetics: Optional[CellLineAlterationProfile] = None
    species_id: str = "human"
    species_scientific_name: str = "Homo sapiens"
    species_common_name: str = "Human"
    human_ortholog_symbol: Optional[str] = None
    nuclease_id: str = "spcas9"
    nuclease_label: str = "SpCas9"
    pam_label: str = "NGG"
    spacer_length_nt: int = 20

    def to_dict(self) -> dict[str, Any]:
        """Serialize the complete report to a JSON-compatible dictionary."""
        return asdict(self)
