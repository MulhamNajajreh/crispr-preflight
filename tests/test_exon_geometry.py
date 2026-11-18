"""Unit tests verifying 5'->3' coordinate math on +1 and -1 strands, 55-nt EJC rule, symmetry, and PAM counting."""

from __future__ import annotations

import unittest

from crispr_preflight.clients.ensembl import (
    GRCH38_NC_MAP,
    EnsemblClient,
    deduplicate_protein_domains,
)
from crispr_preflight.engine.exon_scorer import (
    count_spcas9_cut_sites,
    design_crispri_tss_guides,
    design_exon_sgrnas,
    find_internal_inframe_mets,
    reverse_complement,
    score_transcript_exons,
    select_coding_isoforms,
)
from crispr_preflight.engine.primers import calc_primer_tm, design_tide_ice_primers
from crispr_preflight.models import ProteinDomain


class TestExonGeometryAndRules(unittest.TestCase):
    def test_reverse_strand_ordering_utr_stripping_and_55nt_rule(self) -> None:
        """Verify (-1) reverse-strand 5'->3' ordering, UTR exclusion, symmetry, and final 55-nt NMD escape."""
        # Construct a synthetic 4-exon reverse-strand (-1) transcript:
        # Exon 1 (5' UTR + first 90 bp of CDS): genomic [9000, 9150] (151 bp; CDS starts at 9089 -> 90 bp coding)
        # Exon 2 (internal asymmetric exon, 130 bp, overlaps catalytic domain): genomic [7000, 7129]
        # Exon 3 (internal symmetric penultimate exon, 120 bp -> div3, last 55 bp in EJC-blind zone): genomic [5000, 5119]
        # Exon 4 (last exon, 80 bp coding + 220 bp 3' UTR): genomic [2000, 2299] (CDS ends at 2220 -> [2220, 2299] = 80 bp coding)
        primary_tx = {
            "id": "ENST00000999001",
            "display_name": "TEST1-201",
            "biotype": "protein_coding",
            "strand": -1,
            "is_canonical": 1,
            "MANE": [{"type": "MANE_Select"}],
            "TSL": "tsl1",
            "Translation": {"id": "ENSP00000999001", "start": 2220, "end": 9089, "length": 139},
            "Exon": [
                # Intentionally unsorted to ensure reverse-strand sorting works!
                {"id": "ENSE_E3", "start": 5000, "end": 5119},
                {"id": "ENSE_E1", "start": 9000, "end": 9150},
                {"id": "ENSE_E4", "start": 2000, "end": 2299},
                {"id": "ENSE_E2", "start": 7000, "end": 7129},
            ],
        }
        # Alternative isoform that skips Exon 3
        alt_tx = {
            "id": "ENST00000999002",
            "display_name": "TEST1-202",
            "biotype": "protein_coding",
            "strand": -1,
            "TSL": "tsl1",
            "Translation": {"id": "ENSP00000999002", "start": 2220, "end": 9089, "length": 99},
            "Exon": [
                {"id": "ENSE_E1", "start": 9000, "end": 9150},
                {"id": "ENSE_E2", "start": 7000, "end": 7129},
                {"id": "ENSE_E4", "start": 2000, "end": 2299},
            ],
        }
        domains = [
            ProteinDomain(
                source_db="Pfam",
                accession="PF00069",
                description="Protein kinase domain",
                aa_start=35,
                aa_end=72,
            )
        ]
        # Build a CDS with an internal ATG at codon 18 (within first 150 nt, upstream of domain at aa 35)
        codons = ["GCT"] * 140
        codons[0] = "ATG"
        codons[17] = "ATG"  # Met18!
        cds_seq = "".join(codons)

        audits = score_transcript_exons(
            chrom="17",
            primary_tx=primary_tx,
            alt_coding_txs=[alt_tx],
            domains=domains,
            cds_seq=cds_seq,
        )

        self.assertEqual([a.exon_number for a in audits], [1, 2, 3, 4])
        self.assertEqual([a.exon_id for a in audits], ["ENSE_E1", "ENSE_E2", "ENSE_E3", "ENSE_E4"])

        # Exon 1: coding span [9000, 9089] = 90 bp, triggers 5' Met18 re-initiation warning
        e1 = audits[0]
        self.assertEqual(e1.coding_length_bp, 90)
        self.assertEqual(e1.aa_start, 1)
        self.assertEqual(e1.aa_end, 30)
        self.assertFalse(e1.is_symmetric)  # Exon 1 is first exon with 5' UTR, not internal skippable
        self.assertEqual(e1.reinitiation_met_codon, 18)
        self.assertEqual(e1.nmd_status, "REINIT_5PRIME")

        # Exon 2: internal asymmetric (130 bp % 3 != 0), constitutive, NMD-competent, overlaps Pfam domain -> Score 100!
        e2 = audits[1]
        self.assertEqual(e2.coding_length_bp, 130)
        self.assertFalse(e2.is_symmetric)
        self.assertTrue(e2.is_constitutive)
        self.assertEqual(e2.nmd_status, "COMPETENT")
        self.assertEqual(e2.priority_score, 100)
        self.assertEqual(e2.tier, "RECOMMENDED")
        self.assertEqual(e2.crispick_target_region, "chr17:7000-7129")

        # Exon 3: internal symmetric (120 bp % 3 == 0), skipped in TEST1-202, AND last 55/120 bp (~45.8%) in 55-nt zone
        e3 = audits[2]
        self.assertTrue(e3.is_symmetric)
        self.assertFalse(e3.is_constitutive)
        self.assertIn("TEST1-202", e3.skipped_in_isoforms)
        self.assertAlmostEqual(e3.nmd_escape_fraction, 55 / 120, places=2)
        self.assertEqual(e3.tier, "AVOID")

        # Exon 4: last exon -> 100% NMD escape
        e4 = audits[3]
        self.assertEqual(e4.coding_length_bp, 80)
        self.assertEqual(e4.nmd_status, "ESCAPE_LAST_EXON")
        self.assertEqual(e4.nmd_escape_fraction, 1.0)
        self.assertEqual(e4.tier, "AVOID")

    def test_single_exon_intronless_gene_global_nmd_escape(self) -> None:
        """Intronless / 1-exon transcripts (like SOX2) have 0 splice junctions and must flag ESCAPE_SINGLE_EXON."""
        single_tx = {
            "id": "ENST00000325404",
            "display_name": "SOX2-201",
            "biotype": "protein_coding",
            "strand": 1,
            "is_canonical": 1,
            "Exon": [{"id": "ENSE_SOX2_1", "start": 1000, "end": 3500}],
            "Translation": {"id": "ENSP_SOX2", "start": 1400, "end": 2353, "length": 317},
        }
        domains = [
            ProteinDomain("Pfam", "PF00505", "HMG (high mobility group) box", 40, 110)
        ]
        audits = score_transcript_exons("3", single_tx, [], domains, cds_seq="ATG" + "GCC" * 316)
        self.assertEqual(len(audits), 1)
        self.assertEqual(audits[0].nmd_status, "ESCAPE_SINGLE_EXON")
        self.assertEqual(audits[0].nmd_escape_fraction, 1.0)

    def test_spcas9_pam_cut_site_counting(self) -> None:
        """Verify both sense (NGG) and antisense (CCN) SpCas9 cut sites inside coding interval."""
        # Build a sequence with left_pad=25, coding_len=30, right_pad=25 (total 80 nt of 'A's)
        seq = list("A" * 80)
        # Place a sense NGG at indices 45..47 ('AGG') -> 23-mer starts at i=25, cut is at i+17 = 42 (inside [25, 55))
        seq[45:48] = list("AGG")
        # Place an antisense CCN at indices 24..26 ('CCA') -> 23-mer starts at i=24, cut is at i+6 = 30 (inside [25, 55))
        seq[24:27] = list("CCA")
        # Place an NGG far in the right intron at indices 72..74 -> cut at 69 (outside [25, 55), should NOT count!)
        seq[72:75] = list("TGG")

        cuts = count_spcas9_cut_sites("".join(seq), left_pad=25, coding_len=30)
        self.assertEqual(cuts, 2)

    def test_select_coding_isoforms_keeps_single_exon_primary(self) -> None:
        """Ensure MANE Select / single-exon genes with TSL=None are selected cleanly."""
        gene = {
            "id": "ENSG00000181449",
            "display_name": "SOX2",
            "Transcript": [
                {
                    "id": "ENST00000325404",
                    "display_name": "SOX2-201",
                    "biotype": "protein_coding",
                    "is_canonical": 1,
                    "MANE": [{"type": "MANE_Select"}],
                    "TSL": None,
                    "Exon": [{"id": "E1", "start": 100, "end": 2000}],
                    "Translation": {"id": "ENSP1", "start": 200, "end": 1150, "length": 317},
                }
            ],
        }
        _, summ, alts, _ = select_coding_isoforms(gene)
        self.assertTrue(summ.is_mane_select)
        self.assertEqual(summ.num_exons, 1)
        self.assertEqual(len(alts), 0)

    def test_find_internal_inframe_mets(self) -> None:
        cds = "ATG" + "GCT" * 10 + "ATG" + "GCT" * 5 + "ATG"
        # Codons: 1=ATG, 2..11=GCT, 12=ATG, 13..17=GCT, 18=ATG
        self.assertEqual(find_internal_inframe_mets(cds), [12, 18])

    def test_design_exon_sgrnas_sense_antisense_oligos_and_tttt_filter(self) -> None:
        """Verify 20-nt spacer + PAM extraction, reverse-complement antisense math, TTTT filtering, and U6 oligos."""
        # Build a 90-nt window: left_pad=25, coding_len=40 ([1000..1039]), right_pad=25
        seq = list("A" * 90)
        # 1. Valid '+' strand guide starting at index 20 (cut at 20+17 = 37, inside [25, 65)):
        #    Spacer starts with 'G' (20 nt, 50% GC, ends with 'G') + PAM 'CGG'
        plus_spacer = "GCATGCATGCATGCAATGCG"
        seq[20:43] = list(plus_spacer + "CGG")

        # 2. Valid '-' strand guide starting at index 44 (cut at 44+6 = 50, inside [25, 65)):
        #    Desired 5'->3' spacer on '-' strand starts with 'A' (so CACCG... overhang is tested) + PAM 'AGG'
        minus_spacer = "ACGTACGTACGTACGTCGAG"
        minus_23_plus_strand = reverse_complement(minus_spacer + "AGG")  # starts with 'CCT'
        seq[44:67] = list(minus_23_plus_strand)

        # 3. Invalid '+' strand guide containing Pol III 'TTTT' terminator at index 67 ( wait: cut at 67+17=84 outside; place at index 26? No, don't overwrite!)
        # Let's make window 120 nt with coding_len=70 ([25, 95)) and place TTTT guide at index 68 (cut at 85 inside [25, 95))
        seq.extend(list("A" * 30))
        tttt_spacer = "GCATGCTTTTATGCATGCag".upper()
        seq[68:91] = list(tttt_spacer + "AGG")

        guides = design_exon_sgrnas(
            padded_forward_seq="".join(seq),
            left_pad=25,
            coding_len=70,
            chrom="chr19",
            coding_genomic_start=1000,
            coding_genomic_end=1069,
            gene_strand=-1,  # Reverse-strand gene to test Antisense/Sense & aa_cut math!
            cds_offset_start=31,  # Starts at codon 11
            exon_number=3,
            domains=[ProteinDomain("Pfam", "PF00439", "Bromodomain", 10, 35)],
            gene_symbol="BRD4",
            max_guides=5,
        )

        spacers = {g.spacer_20nt: g for g in guides}
        self.assertIn(plus_spacer, spacers)
        self.assertIn(minus_spacer, spacers)
        self.assertNotIn(tttt_spacer, spacers)  # TTTT Pol III terminator strictly excluded!

        g_plus = spacers[plus_spacer]
        self.assertEqual(g_plus.pam, "CGG")
        self.assertEqual(g_plus.genomic_strand, "+")
        self.assertEqual(g_plus.transcript_orientation, "Antisense")  # + genomic on -1 gene
        self.assertEqual(g_plus.cloning_fwd_oligo, f"CACC{plus_spacer}")
        self.assertEqual(g_plus.cloning_rev_oligo, f"AAAC{reverse_complement(plus_spacer)}")
        self.assertEqual(g_plus.overlaps_domain, "Bromodomain")

        g_minus = spacers[minus_spacer]
        self.assertEqual(g_minus.pam, "AGG")
        self.assertEqual(g_minus.genomic_strand, "-")
        self.assertEqual(g_minus.transcript_orientation, "Sense")  # - genomic on -1 gene
        self.assertEqual(g_minus.cloning_fwd_oligo, f"CACCG{minus_spacer}")
        self.assertEqual(g_minus.cloning_rev_oligo, f"AAAC{reverse_complement(minus_spacer)}C")

    def test_deduplicate_protein_domains_collapses_overlapping_databases(self) -> None:
        """Overlapping SuperFamily/Smart/CDD hits for the same Pfam fold must collapse cleanly while keeping tandem repeats."""
        raw_domains = [
            ProteinDomain("SuperFamily", "SSF47370", "Bromodomain-like", 55, 168),
            ProteinDomain("Smart", "SM00297", "Bromo", 62, 160),
            ProteinDomain("Pfam", "PF00439", "Bromodomain", 75, 147),
            ProteinDomain("Pfam", "PF00439", "Bromodomain", 368, 440),
            ProteinDomain("CDD", "cd05510", "BrD_BRD4_I", 352, 456),
            ProteinDomain("Pfam", "PF17035", "BET extra-terminal domain", 601, 683),
        ]
        deduped = deduplicate_protein_domains(raw_domains)
        self.assertEqual(len(deduped), 3)
        self.assertEqual([d.source_db for d in deduped], ["Pfam", "Pfam", "Pfam"])
        self.assertEqual([(d.aa_start, d.aa_end) for d in deduped], [(75, 147), (368, 440), (601, 683)])

    def test_ebi_coordinates_fallback_reverse_strand_normalization(self) -> None:
        """Verify _build_gene_from_ebi_coords normalizes reverse-strand begin > end into start <= end."""
        gc_list = [
            {
                "ensemblGeneId": "ENSG00000142208",
                "ensemblTranscriptId": "ENST00000554581.5",
                "ensemblTranslationId": "ENSP00000451828.1",
                "genomicLocation": {
                    "assemblyName": "GRCh38",
                    "chromosome": "14",
                    "reverseStrand": True,
                    "exon": [
                        {
                            "id": "ENSE001",
                            "genomeLocation": {"begin": {"position": 104792662}, "end": {"position": 104792538}},
                            "proteinLocation": {"begin": {"position": 1}, "end": {"position": 42}},
                        },
                        {
                            "id": "ENSE002",
                            "genomeLocation": {"begin": {"position": 104780216}, "end": {"position": 104780088}},
                            "proteinLocation": {"begin": {"position": 43}, "end": {"position": 480}},
                        },
                    ],
                },
            }
        ]
        gene_obj = EnsemblClient._build_gene_from_ebi_coords(
            ensembl_id="ENSG00000142208",
            gene_symbol="AKT1",
            approved_name="AKT serine/threonine kinase 1",
            gc_list=gc_list,
            canon_seq_len=480,
        )
        self.assertIsNotNone(gene_obj)
        assert gene_obj is not None
        self.assertEqual(gene_obj["strand"], -1)
        self.assertEqual(gene_obj["seq_region_name"], "14")
        tx = gene_obj["Transcript"][0]
        self.assertEqual(tx["id"], "ENST00000554581")
        self.assertEqual(tx["Translation"]["length"], 480)
        for ex in tx["Exon"]:
            self.assertLessEqual(ex["start"], ex["end"])
        self.assertEqual(len(GRCH38_NC_MAP), 25)

    def test_santalucia_tm_and_tide_ice_primer_designer(self) -> None:
        """Verify SantaLucia nearest-neighbor Tm and TIDE/ICE PCR amplicon geometry (430-740 bp)."""
        tm = calc_primer_tm("GCTAGCTAGCTAGCTAGCTA")  # 20-nt, 50% GC in standard PCR buffer
        self.assertGreater(tm, 55.0)
        self.assertLess(tm, 65.0)

        # Construct a realistic 900-bp mixed-GC flanking window [7675000..7675899] with cut at 7675450
        motif = "ACGTGCATGCTAGCTGATCGTAGCTAGCATCG"  # 32 bp, 50% GC, balanced, non-palindromic 3' steps
        flank = (motif * 29)[:900]
        pair = design_tide_ice_primers(
            window_plus_seq=flank,
            window_genomic_start=7675000,
            cut_genomic_pos=7675450,
            chrom="chr17",
        )
        self.assertIsNotNone(pair)
        assert pair is not None
        self.assertGreaterEqual(pair.amplicon_bp, 430)
        self.assertLessEqual(pair.amplicon_bp, 740)
        self.assertGreaterEqual(pair.fwd_tm, 57.0)
        self.assertLessEqual(pair.fwd_tm, 63.5)
        self.assertLessEqual(abs(pair.fwd_tm - pair.rev_tm), 2.8)
        self.assertIn(pair.recommended_seq_primer, ("Forward", "Reverse"))

    def test_design_crispri_tss_guides_forward_and_reverse_strand(self) -> None:
        """Verify Horlbeck/Weissman [-50, +300 bp] promoter CRISPRi window on +1 and -1 strands."""
        # Build 900-bp promoter slice [10000..10899] with TSS at 10450
        seq = list("A" * 900)
        sp1 = "GCATGCATGCATGCAATGCG"
        seq[500:523] = list(sp1 + "AGG")  # cut at 10517 -> +67 bp downstream of +1 strand TSS (10450)

        guides_plus = design_crispri_tss_guides(
            tss_window_plus_seq="".join(seq),
            window_genomic_start=10000,
            tss_pos=10450,
            chrom="chr19",
            gene_strand=1,
            gene_symbol="PLK1",
        )
        self.assertTrue(any(g.spacer_20nt == sp1 for g in guides_plus))
        g1 = next(g for g in guides_plus if g.spacer_20nt == sp1)
        self.assertEqual(g1.mode, "CRISPRi")
        self.assertEqual(g1.exon_number, 0)
        self.assertGreaterEqual(g1.tss_offset_bp or -999, -50)
        self.assertLessEqual(g1.tss_offset_bp or 999, 300)

        # Reverse-strand (-1) gene with TSS at 10550: downstream (+60 bp) is at LOWER genomic coordinate (10490 -> index 490)
        seq_rev = list("A" * 900)
        seq_rev[473:496] = list(sp1 + "CGG")  # cut at 10000 + 490 = 10490 -> offset = 10550 - 10490 = +60 bp!
        guides_minus = design_crispri_tss_guides(
            tss_window_plus_seq="".join(seq_rev),
            window_genomic_start=10000,
            tss_pos=10550,
            chrom="chr19",
            gene_strand=-1,
            gene_symbol="BRD4",
        )
        self.assertTrue(any(g.spacer_20nt == sp1 for g in guides_minus))
        gm = next(g for g in guides_minus if g.spacer_20nt == sp1)
        self.assertEqual(gm.tss_offset_bp, 60)


if __name__ == "__main__":
    unittest.main()
