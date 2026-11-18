"""Unit tests verifying 5'->3' coordinate math on +1 and -1 strands, 55-nt EJC rule, symmetry, and PAM counting."""

from __future__ import annotations

import unittest

from crispr_preflight.clients.ensembl import (
    GRCH38_NC_MAP,
    EnsemblClient,
    deduplicate_protein_domains,
)
from crispr_preflight.clients.multispecies import MultiSpeciesClient
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
from crispr_preflight.nucleases import NUCLEASE_REGISTRY, resolve_nuclease
from crispr_preflight.species import SPECIES_REGISTRY, resolve_species


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

    def test_multispecies_registry_and_mygene_zero_based_conversion(self) -> None:
        """Verify all 6 organism aliases resolve cleanly and 0-based half-open MyGene exons convert to 1-based closed."""
        self.assertEqual(set(SPECIES_REGISTRY.keys()), {"human", "mouse", "rat", "zebrafish", "fly", "worm"})
        self.assertEqual(resolve_species("Mus musculus").assembly, "GRCm39")
        self.assertEqual(resolve_species("danio_rerio").assembly, "GRCz11")
        self.assertEqual(resolve_species("drosophila").assembly, "BDGP6.46")
        self.assertEqual(resolve_species("c. elegans").assembly, "WBcel235")
        self.assertEqual(resolve_species("rn").assembly, "GRCr8")
        with self.assertRaises(ValueError):
            resolve_species("unicorn")

        # Synthetic 2-exon mouse MyGene hit (0-based half-open coordinates):
        # Exon 1: [1000, 1200) -> 1-based closed [1001, 1200] (200 bp); CDS starts at 1050 -> 1-based 1051 (150 bp coding)
        # Exon 2: [2000, 2300) -> 1-based closed [2001, 2300] (300 bp); CDS ends at 2150 -> 1-based 2150 (150 bp coding)
        # Total coding = 300 bp -> 99 aa + stop codon!
        hit = {
            "symbol": "Trp53",
            "name": "transformation related protein 53",
            "ensembl": [
                {"gene": "ENSMSIG00000012345"},  # Non-reference mouse strain should be skipped in favor of ENSMUSG!
                {"gene": "ENSMUSG00000059552"},
            ],
            "exons": {
                "NM_011640": {
                    "chr": "11",
                    "strand": 1,
                    "txstart": 1000,
                    "txend": 2300,
                    "cdsstart": 1050,
                    "cdsend": 2150,
                    "position": [[1000, 1200], [2000, 2300]],
                }
            },
        }
        mouse_cfg = resolve_species("mouse")
        ens_id = MultiSpeciesClient._extract_ensembl_gene_id(hit["ensembl"], "mouse")
        self.assertEqual(ens_id, "ENSMUSG00000059552")

        gene_obj, tss_pos = MultiSpeciesClient.build_gene_from_mygene_exons(
            hit=hit,
            gene_id=ens_id or "ENSMUSG00000059552",
            gene_symbol="Trp53",
            description="transformation related protein 53",
            species_cfg=mouse_cfg,
        )
        self.assertIsNotNone(gene_obj)
        assert gene_obj is not None
        self.assertEqual(gene_obj["id"], "ENSMUSG00000059552")
        self.assertEqual(gene_obj["assembly_name"], "GRCm39")
        self.assertEqual(tss_pos, 1001)
        tx = gene_obj["Transcript"][0]
        self.assertEqual(tx["Translation"]["start"], 1051)
        self.assertEqual(tx["Translation"]["end"], 2150)
        self.assertEqual(tx["Translation"]["length"], 99)
        self.assertEqual([(e["start"], e["end"]) for e in tx["Exon"]], [(1001, 1200), (2001, 2300)])

    def test_multi_nuclease_orthologs_spg_hifi_sacas9_and_cas12a(self) -> None:
        """Verify SpG (NGN), SpCas9-HF (native 5'-G), SaCas9 (21-nt + NNGRRT), and Cas12a (5'-TTTV + 23-nt)."""
        self.assertEqual(set(NUCLEASE_REGISTRY.keys()), {"spcas9", "spcas9_hf", "spg", "sacas9", "cas12a"})
        self.assertEqual(resolve_nuclease("aav").nuclease_id, "sacas9")
        self.assertEqual(resolve_nuclease("cpf1").nuclease_id, "cas12a")
        self.assertEqual(resolve_nuclease("hf1").nuclease_id, "spcas9_hf")
        self.assertEqual(resolve_nuclease("ngn").nuclease_id, "spg")
        with self.assertRaises(ValueError):
            resolve_nuclease("meganuclease_xyz")

        # 1. SpG (20-nt + NGN) rescues exons that have NGA/NGC/NGT but zero NGG!
        seq_spg = list("A" * 100)
        seq_spg[30:53] = list(("A" * 20) + "AGA")  # Single G in entire 100-bp window -> PAM = AGA (1 NGN site, 0 NGG!)
        self.assertEqual(count_spcas9_cut_sites("".join(seq_spg), 25, 50, nuclease_id="spcas9"), 0)
        self.assertEqual(count_spcas9_cut_sites("".join(seq_spg), 25, 50, nuclease_id="spg"), 1)

        # 2. SpCas9-HF ranks native 5'-G spacer strictly above non-G spacer of identical GC%
        seq_hf = list("A" * 140)
        g_start_sp = "GCATGCATGCATGCAATGCG"  # starts with G (50% GC)
        c_start_sp = "CATGCATGCATGCAATGCGC"  # starts with C (50% GC)
        seq_hf[25:48] = list(g_start_sp + "CGG")
        seq_hf[65:88] = list(c_start_sp + "CGG")
        hf_guides = design_exon_sgrnas(
            padded_forward_seq="".join(seq_hf),
            left_pad=20,
            coding_len=90,
            chrom="chr17",
            coding_genomic_start=1000,
            coding_genomic_end=1089,
            gene_strand=1,
            cds_offset_start=1,
            exon_number=2,
            domains=[],
            gene_symbol="TP53",
            nuclease_id="spcas9_hf",
        )
        self.assertGreaterEqual(len(hf_guides), 2)
        self.assertEqual(hf_guides[0].spacer_20nt, g_start_sp)
        self.assertGreater(hf_guides[0].sequence_score, hf_guides[1].sequence_score + 15)

        # 3. SaCas9 (21-nt + 3'-NNGRRT) on '+' and '-' strands
        seq_sa = list("A" * 140)
        sa_plus = "GCATGCATGCATGCAATGCGA"  # 21-nt
        sa_minus = "ACGTACGTACGTACGTCGAGA"  # 21-nt
        seq_sa[22:49] = list(sa_plus + "CTGGGT")  # NNGRRT = CTGGGT (cut at 22+18=40 inside [25, 115))
        seq_sa[65:92] = list(reverse_complement(sa_minus + "TCGAAT"))  # NNGRRT = TCGAAT (cut at 65+9=74)
        sa_guides = design_exon_sgrnas(
            padded_forward_seq="".join(seq_sa),
            left_pad=25,
            coding_len=90,
            chrom="chr11",
            coding_genomic_start=2000,
            coding_genomic_end=2089,
            gene_strand=1,
            cds_offset_start=1,
            exon_number=4,
            domains=[],
            gene_symbol="Trp53",
            nuclease_id="sacas9",
        )
        sa_map = {g.spacer_20nt: g for g in sa_guides}
        self.assertIn(sa_plus, sa_map)
        self.assertIn(sa_minus, sa_map)
        self.assertEqual(sa_map[sa_plus].spacer_length, 21)
        self.assertEqual(sa_map[sa_plus].pam, "CTGGGT")
        self.assertEqual(sa_map[sa_minus].pam, "TCGAAT")

        # 4. Cas12a (5'-TTTV + 23-nt) on '+' and '-' strands + AGAT/AAAA BsmBI crRNA oligos
        seq_c12 = list("A" * 140)
        c12_plus = "CGATGCATGCATGCAATGCGATC"  # 23-nt
        c12_minus = "GCGTACGTACGTACGTCGAGATC"  # 23-nt
        seq_c12[20:47] = list("TTTC" + c12_plus)  # cut at 20+23 = 43 inside [25, 115)
        seq_c12[65:92] = list(reverse_complement("TTTG" + c12_minus))  # cut at 65+4 = 69 inside [25, 115)
        c12_guides = design_exon_sgrnas(
            padded_forward_seq="".join(seq_c12),
            left_pad=25,
            coding_len=90,
            chrom="chrI",
            coding_genomic_start=3000,
            coding_genomic_end=3089,
            gene_strand=1,
            cds_offset_start=1,
            exon_number=1,
            domains=[],
            gene_symbol="cep-1",
            nuclease_id="cas12a",
        )
        c12_map = {g.spacer_20nt: g for g in c12_guides}
        self.assertIn(c12_plus, c12_map)
        self.assertIn(c12_minus, c12_map)
        g_c12 = c12_map[c12_plus]
        self.assertEqual(g_c12.pam, "TTTC")
        self.assertEqual(g_c12.pam_orientation, "5prime")
        self.assertEqual(g_c12.spacer_length, 23)
        self.assertEqual(g_c12.cloning_fwd_oligo, f"AGAT{c12_plus}")
        self.assertEqual(g_c12.cloning_rev_oligo, f"AAAA{reverse_complement(c12_plus)}")


if __name__ == "__main__":
    unittest.main()
