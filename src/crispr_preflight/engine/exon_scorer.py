"""Exon-by-exon coordinate geometry, isoform constitutiveness, NMD escape, symmetry, and PAM scoring.

Literature basis:
- Smits et al., Nature Methods (2019): Residual protein expression in frameshift KO lines.
- Tuladhar et al., Nature Communications (2019): Symmetric exon skipping & downstream AUG re-initiation.
- Lindeboom et al., Nature Genetics (2019) [NMDetective]: 55-nt EJC rule, start-proximal & >400-nt exon rules.
- Shi et al., Nature Biotechnology (2015): Domain-focused CRISPR mutagenesis defeats in-frame (-3 bp) escapers.
"""

from __future__ import annotations

import re
from typing import Any, Optional

from crispr_preflight.engine.primers import design_tide_ice_primers
from crispr_preflight.models import ExonAudit, GuideCandidate, IsoformSummary, ProteinDomain
from crispr_preflight.nucleases import NucleaseConfig, resolve_nuclease

_DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(seq: str) -> str:
    """Return exact 5'->3' reverse complement of an uppercase A/C/G/T DNA string."""
    return seq.upper().translate(_DNA_COMPLEMENT)[::-1]


def parse_tsl(raw_tsl: Any) -> Optional[str]:
    """Normalize Ensembl TSL field (e.g., 'tsl1 (assigned...)' or {'value': 1}) to '1'..'5' or None."""
    if raw_tsl is None:
        return None
    if isinstance(raw_tsl, dict):
        val = raw_tsl.get("value")
        return str(val) if val is not None else None
    text = str(raw_tsl).strip().lower()
    m = re.match(r"^(?:tsl)?([1-5])\b", text)
    return m.group(1) if m else None


def extract_mane_types(transcript: dict[str, Any]) -> set[str]:
    """Extract MANE tags ('MANE_Select', 'MANE_Plus_Clinical') from an Ensembl transcript dict."""
    tags: set[str] = set()
    mane_field = transcript.get("MANE")
    if isinstance(mane_field, list):
        for item in mane_field:
            if isinstance(item, dict) and item.get("type"):
                tags.add(str(item["type"]))
    elif isinstance(mane_field, dict) and mane_field.get("type"):
        tags.add(str(mane_field["type"]))
    if transcript.get("is_mane_select"):
        tags.add("MANE_Select")
    return tags


def select_coding_isoforms(
    gene_data: dict[str, Any],
) -> tuple[dict[str, Any], IsoformSummary, list[dict[str, Any]], list[IsoformSummary]]:
    """Select the primary coding transcript (MANE Select > Ensembl Canonical > longest) and well-supported alternatives."""
    raw_transcripts = gene_data.get("Transcript") or []
    coding_entries: list[tuple[dict[str, Any], IsoformSummary]] = []

    for tx in raw_transcripts:
        if tx.get("biotype") != "protein_coding":
            continue
        trans = tx.get("Translation")
        if not isinstance(trans, dict) or trans.get("start") is None or trans.get("end") is None:
            continue
        mane_tags = extract_mane_types(tx)
        is_mane = "MANE_Select" in mane_tags
        is_mane_plus = "MANE_Plus_Clinical" in mane_tags
        is_canon = bool(tx.get("is_canonical"))
        tsl = parse_tsl(tx.get("TSL") or tx.get("tsl"))
        prot_len = int(trans.get("length") or 0)
        exons = tx.get("Exon") or []
        cds_len = _compute_transcript_cds_bp(exons, int(trans["start"]), int(trans["end"]))

        summary = IsoformSummary(
            transcript_id=str(tx["id"]),
            display_name=str(tx.get("display_name") or tx["id"]),
            is_mane_select=is_mane,
            is_mane_plus_clinical=is_mane_plus,
            is_canonical=is_canon,
            tsl=tsl,
            cds_length_nt=cds_len,
            protein_length_aa=prot_len,
            num_exons=len(exons),
            protein_id=str(trans.get("id")) if trans.get("id") else None,
        )
        coding_entries.append((tx, summary))

    if not coding_entries:
        symbol = gene_data.get("display_name") or gene_data.get("id")
        raise ValueError(
            f"Gene '{symbol}' has no annotated protein_coding transcripts with a CDS in Ensembl GRCh38."
        )

    # Sort: MANE Select first, then Canonical, then TSL 1, then longest protein
    coding_entries.sort(
        key=lambda pair: (
            0 if pair[1].is_mane_select else 1,
            0 if pair[1].is_canonical else 1,
            0 if pair[1].tsl == "1" else 1,
            -pair[1].protein_length_aa,
        )
    )

    primary_raw, primary_summary = coding_entries[0]

    # Keep high-confidence alternative isoforms for constitutive exon comparison:
    # In modern Ensembl/GENCODE (v46+), `gencode_primary == 1` + MANE marks biologically
    # expressed primary/alternative isoforms (avoiding dozens of partial computational fragments).
    alt_raws: list[dict[str, Any]] = []
    alt_summaries: list[IsoformSummary] = []
    has_gencode_primary_flag = any("gencode_primary" in tx for tx, _ in coding_entries)
    min_credible_aa = max(50, int(primary_summary.protein_length_aa * 0.35))

    for tx, summ in coding_entries[1:]:
        if has_gencode_primary_flag:
            is_well_supported = (
                bool(tx.get("gencode_primary"))
                or summ.is_mane_select
                or summ.is_mane_plus_clinical
                or summ.is_canonical
            )
        else:
            is_well_supported = (
                summ.is_mane_plus_clinical
                or summ.is_canonical
                or summ.tsl in ("1", "2")
                or (summ.tsl is None and summ.protein_length_aa >= min_credible_aa)
            )
        if is_well_supported:
            alt_raws.append(tx)
            alt_summaries.append(summ)

    return primary_raw, primary_summary, alt_raws, alt_summaries


def _compute_transcript_cds_bp(exons: list[dict[str, Any]], t_pos1: int, t_pos2: int) -> int:
    t_start, t_end = min(t_pos1, t_pos2), max(t_pos1, t_pos2)
    total = 0
    for ex in exons:
        s, e = min(int(ex["start"]), int(ex["end"])), max(int(ex["start"]), int(ex["end"]))
        ov_s, ov_e = max(s, t_start), min(e, t_end)
        if ov_s <= ov_e:
            total += ov_e - ov_s + 1
    return total


def find_internal_inframe_mets(cds_seq: str, max_nt: int = 450) -> list[int]:
    """Return 1-based amino-acid positions of internal in-frame ATG (Met) codons within first max_nt."""
    if not cds_seq or len(cds_seq) < 6:
        return []
    seq = cds_seq.upper()
    limit_codons = min(len(seq) // 3, max_nt // 3)
    mets: list[int] = []
    for codon_num in range(2, limit_codons + 1):
        codon = seq[(codon_num - 1) * 3 : codon_num * 3]
        if codon == "ATG":
            mets.append(codon_num)
    return mets


def _scan_nuclease_sites(
    seq: str,
    p_start: int,
    c_start: int,
    c_end: int,
    nuc: NucleaseConfig,
) -> list[tuple[str, str, str, int, int, int]]:
    """Scan a forward genomic slice for valid (spacer, pam, strand, cut_pos, g_win_start, g_win_end) sites."""
    valid_bases = {"A", "C", "G", "T"}
    nid = nuc.nuclease_id
    w_len = nuc.window_length
    n = len(seq)
    hits: list[tuple[str, str, str, int, int, int]] = []
    if n < w_len:
        return hits

    for i in range(0, n - w_len + 1):
        window = seq[i : i + w_len]
        if any(ch not in valid_bases for ch in window):
            continue
        g_ws = p_start + i
        g_we = g_ws + w_len - 1

        if nid in ("spcas9", "spcas9_hf"):
            # 20-nt + 3'-NGG (23 bp)
            if window[21:23] == "GG":
                cut_idx = i + 17
                if c_start <= cut_idx < c_end:
                    hits.append((window[0:20], window[20:23], "+", p_start + cut_idx, g_ws, g_we))
            if window[0:2] == "CC":
                cut_idx = i + 6
                if c_start <= cut_idx < c_end:
                    rc_w = reverse_complement(window)
                    hits.append((rc_w[0:20], rc_w[20:23], "-", p_start + cut_idx, g_ws, g_we))

        elif nid == "spg":
            # 20-nt + 3'-NGN (23 bp; Walton et al., Science 2020)
            if window[21] == "G":
                cut_idx = i + 17
                if c_start <= cut_idx < c_end:
                    hits.append((window[0:20], window[20:23], "+", p_start + cut_idx, g_ws, g_we))
            if window[1] == "C":
                cut_idx = i + 6
                if c_start <= cut_idx < c_end:
                    rc_w = reverse_complement(window)
                    hits.append((rc_w[0:20], rc_w[20:23], "-", p_start + cut_idx, g_ws, g_we))

        elif nid == "sacas9":
            # 21-nt + 3'-NNGRRT (27 bp; Ran et al., Nature 2015)
            if window[23] == "G" and window[24] in ("A", "G") and window[25] in ("A", "G") and window[26] == "T":
                cut_idx = i + 18
                if c_start <= cut_idx < c_end:
                    hits.append((window[0:21], window[21:27], "+", p_start + cut_idx, g_ws, g_we))
            if window[0] == "A" and window[1] in ("C", "T") and window[2] in ("C", "T") and window[3] == "C":
                cut_idx = i + 9
                if c_start <= cut_idx < c_end:
                    rc_w = reverse_complement(window)
                    hits.append((rc_w[0:21], rc_w[21:27], "-", p_start + cut_idx, g_ws, g_we))

        elif nid == "cas12a":
            # 5'-TTTV + 23-nt (27 bp; Zetsche et al., Cell 2015; staggered cut centered at +19 nt of spacer)
            if window[0:3] == "TTT" and window[3] in ("A", "C", "G"):
                cut_idx = i + 23
                if c_start <= cut_idx < c_end:
                    hits.append((window[4:27], window[0:4], "+", p_start + cut_idx, g_ws, g_we))
            if window[24:27] == "AAA" and window[23] in ("C", "G", "T"):
                cut_idx = i + 4
                if c_start <= cut_idx < c_end:
                    rc_w = reverse_complement(window)
                    hits.append((rc_w[4:27], rc_w[0:4], "-", p_start + cut_idx, g_ws, g_we))

    return hits


def count_spcas9_cut_sites(
    padded_forward_seq: str,
    left_pad: int,
    coding_len: int,
    nuclease_id: str = "spcas9",
) -> int:
    """Count nuclease cut sites landing inside [left_pad, left_pad + coding_len)."""
    if not padded_forward_seq or coding_len <= 0:
        return 0
    nuc = resolve_nuclease(nuclease_id)
    hits = _scan_nuclease_sites(
        seq=padded_forward_seq.upper(),
        p_start=1,
        c_start=left_pad,
        c_end=left_pad + coding_len,
        nuc=nuc,
    )
    return len(hits)


def _build_cloning_oligos(spacer: str, nuc: NucleaseConfig) -> tuple[str, str]:
    """Build standard Golden Gate U6 cloning oligos matched to the nuclease backbone."""
    sp = spacer.upper()
    rc = reverse_complement(sp)
    if nuc.nuclease_id == "cas12a":
        # Zhang-lab AsCas12a/LbCas12a U6 direct-repeat BsmBI overhangs (pY010 / pY016 / pRDA_052)
        return f"AGAT{sp}", f"AAAA{rc}"
    if sp.startswith("G"):
        return f"CACC{sp}", f"AAAC{rc}"
    return f"CACCG{sp}", f"AAAC{rc}C"


def _score_nuclease_spacer(
    spacer: str,
    pam: str,
    nuc: NucleaseConfig,
    base_score: int = 82,
) -> tuple[int, int, list[str], str, str]:
    """Evaluate GC%, PAM/seed rules, homopolymers, and Golden Gate restriction sites for any nuclease."""
    s_len = len(spacer)
    gc_count = spacer.count("G") + spacer.count("C")
    gc_pct = int(round((gc_count / float(s_len)) * 100))
    s_score = base_score
    notes: list[str] = []
    nid = nuc.nuclease_id

    opt_low = 35 if nid == "cas12a" else 40
    if opt_low <= gc_pct <= 65:
        s_score += 6
        notes.append(f"Optimal GC ({gc_pct}%)")
    elif 30 <= gc_pct <= 70:
        s_score -= 4
        notes.append(f"Moderate GC ({gc_pct}%)")
    else:
        s_score -= 22
        notes.append(f"Extreme GC ({gc_pct}%)")

    if nuc.pam_orientation == "5prime":
        # Cas12a (5'-TTTV + 23-nt): 5' seed is spacer[0:8] immediately 3' of TTTV PAM (Kim et al., Nat Methods 2017)
        if spacer[0] == "T":
            s_score -= 14
            notes.append("Disfavored T at pos 1 after TTTV")
        elif spacer[0] in ("C", "G"):
            s_score += 5
            notes.append(f"Favored {spacer[0]} at pos 1")
        seed_gc = spacer[0:8].count("G") + spacer[0:8].count("C")
        if 3 <= seed_gc <= 6:
            s_score += 4
            notes.append("Balanced 5' seed")
        if len(pam) >= 4 and pam[3] in ("C", "G"):
            s_score += 2
    else:
        # 3'-PAM Cas9 orthologs (SpCas9, SpCas9-HF, SpG, SaCas9): seed is 3' PAM-proximal end
        last_nt = spacer[-1]
        if last_nt == "G":
            s_score += 6
            notes.append(f"G at pos {s_len}")
        elif last_nt == "C":
            s_score += 2
        elif last_nt == "T":
            s_score -= 8
            notes.append(f"Disfavored T at pos {s_len}")

        if nid in ("spcas9", "spcas9_hf"):
            if pam[0] in ("C", "G"):
                s_score += 2
            elif pam[0] == "T":
                s_score -= 3
        elif nid == "spg":
            # Walton et al. (2020): NGG/NGC/NGA slightly outperform NGT
            if len(pam) >= 3 and pam[2] == "G":
                s_score += 4
                notes.append(f"Canonical {pam} PAM")
            elif len(pam) >= 3 and pam[2] in ("A", "C"):
                s_score += 2
                notes.append(f"Active {pam} PAM")
            elif len(pam) >= 3 and pam[2] == "T":
                s_score -= 4
                notes.append(f"Moderate {pam} PAM")
        elif nid == "sacas9":
            # Ran et al. (2015): NNGGGT / NNGAAT / NNGAGT / NNGGAT
            if len(pam) >= 6 and pam[3:5] == "GG":
                s_score += 4
                notes.append(f"Optimal {pam} PAM")
            elif len(pam) >= 6:
                s_score += 2

        if nid == "spcas9_hf":
            # Kim et al. (2017) & Kulcsar et al. (2017): unmatched 5'-G impairs eSpCas9(1.1) & SpCas9-HF1 under U6
            if spacer.startswith("G"):
                s_score += 12
                notes.append("Native 5'-G (ideal U6 match for eSpCas9/HF1)")
            else:
                s_score -= 14
                notes.append(f"Starts w/ {spacer[0]} (use RNP or tRNA-sgRNA for eSpCas9/HF1)")
        elif nid == "spg" and spacer.startswith("G"):
            s_score += 4
            notes.append("Native 5'-G")

        prox4 = spacer[-4:]
        if "TT" in prox4:
            s_score -= 15
            notes.append("Graf PAM-proximal TT motif")
        if "GCC" in prox4:
            s_score -= 12
            notes.append("Graf PAM-proximal GCC motif")

    if "GGGG" in spacer:
        s_score -= 18
        notes.append("GGGG G-quadruplex risk")
    elif "AAAA" in spacer or "CCCC" in spacer:
        s_score -= 8
        notes.append("4-nt homopolymer")

    fwd_oligo, rev_oligo = _build_cloning_oligos(spacer, nuc)
    if nid == "sacas9":
        if "GGTCTC" in fwd_oligo or "GAGACC" in fwd_oligo:
            s_score -= 6
            notes.append("Internal BsaI site (pX601)")
    elif nid == "cas12a":
        if "CGTCTC" in fwd_oligo or "GAGACG" in fwd_oligo:
            s_score -= 6
            notes.append("Internal BsmBI site")
    else:
        if "CGTCTC" in fwd_oligo or "GAGACG" in fwd_oligo:
            s_score -= 6
            notes.append("Internal BsmBI site")
        if "GAAGAC" in fwd_oligo or "GTCTTC" in fwd_oligo:
            s_score -= 6
            notes.append("Internal BbsI site")

    return gc_pct, s_score, notes, fwd_oligo, rev_oligo


def design_exon_sgrnas(
    padded_forward_seq: str,
    left_pad: int,
    coding_len: int,
    chrom: str,
    coding_genomic_start: int,
    coding_genomic_end: int,
    gene_strand: int,
    cds_offset_start: int,
    exon_number: int,
    domains: list[ProteinDomain],
    gene_symbol: str = "GENE",
    is_single_exon: bool = False,
    max_guides: int = 5,
    nuclease_id: str = "spcas9",
) -> list[GuideCandidate]:
    """Extract and rank clean on-target guide candidates inside a coding exon for the chosen nuclease."""
    if not padded_forward_seq or coding_len <= 0:
        return []

    nuc = resolve_nuclease(nuclease_id)
    seq = padded_forward_seq.upper()
    c_start = left_pad
    c_end = left_pad + coding_len
    p_start = coding_genomic_start - left_pad

    raw_hits = _scan_nuclease_sites(seq, p_start, c_start, c_end, nuc)
    scored: list[GuideCandidate] = []
    seen_spacers: set[tuple[str, int]] = set()

    for spacer, pam, g_strand, cut_pos, g_s_win, g_e_win in raw_hits:
        if "TTTT" in spacer:
            # Hard-exclude Pol III U6 transcription termination signal
            continue
        if (spacer, cut_pos) in seen_spacers:
            continue
        seen_spacers.add((spacer, cut_pos))

        rel_nt = (
            (cut_pos - coding_genomic_start)
            if gene_strand >= 0
            else (coding_genomic_end - cut_pos)
        )
        cds_nt = cds_offset_start + max(0, rel_nt)
        aa_cut = ((cds_nt - 1) // 3) + 1

        dom_hit: Optional[str] = None
        for dom in domains:
            if dom.aa_start <= aa_cut <= dom.aa_end:
                dom_hit = dom.description
                break

        gc_pct, s_score, notes, fwd_oligo, rev_oligo = _score_nuclease_spacer(
            spacer=spacer,
            pam=pam,
            nuc=nuc,
            base_score=82,
        )

        # Conserved domain cut bonus (extra weight for intronless 1-exon genes where NMD is absent)
        if dom_hit:
            s_score += 18 if is_single_exon else 7
            notes.append(f"Cuts in domain ({dom_hit[:24]})")
        elif is_single_exon:
            s_score -= 12
            notes.append("Outside annotated domain on 1-exon gene")

        # Interior vs splice-edge cut
        if coding_len > 10 and (rel_nt < 3 or rel_nt >= coding_len - 3):
            s_score -= 8
            notes.append("Near splice junction (<=3 bp)")
        else:
            s_score += 2

        s_score = max(0, min(100, s_score))
        tx_orient = (
            "Sense"
            if ((g_strand == "+" and gene_strand >= 0) or (g_strand == "-" and gene_strand < 0))
            else "Antisense"
        )

        scored.append(
            GuideCandidate(
                guide_id="",
                exon_number=exon_number,
                spacer_20nt=spacer,
                pam=pam,
                genomic_strand=g_strand,
                transcript_orientation=tx_orient,
                chrom=chrom,
                genomic_cut_pos=cut_pos,
                aa_cut=aa_cut,
                gc_percent=gc_pct,
                overlaps_domain=dom_hit,
                sequence_score=s_score,
                quality_notes=notes,
                cloning_fwd_oligo=fwd_oligo,
                cloning_rev_oligo=rev_oligo,
                mode="KO",
                genomic_start_23bp=g_s_win,
                genomic_end_23bp=g_e_win,
                nuclease=nuc.short_tag,
                spacer_length=nuc.spacer_length,
                pam_orientation=nuc.pam_orientation,
                cloning_backbone_hint=nuc.cloning_backbone,
            )
        )

    scored.sort(key=lambda g: (-g.sequence_score, g.aa_cut, g.genomic_cut_pos))

    # First pass: select spatially non-overlapping cuts (>=5 bp apart) for independent guide diversity
    selected: list[GuideCandidate] = []
    for cand in scored:
        if len(selected) >= max_guides:
            break
        if all(abs(cand.genomic_cut_pos - prev.genomic_cut_pos) >= 5 for prev in selected):
            selected.append(cand)

    # Backfill if exon is very short and had fewer than max_guides spaced sites
    if len(selected) < max_guides:
        chosen_keys = {(g.spacer_20nt, g.genomic_cut_pos) for g in selected}
        for cand in scored:
            if len(selected) >= max_guides:
                break
            key = (cand.spacer_20nt, cand.genomic_cut_pos)
            if key not in chosen_keys:
                selected.append(cand)
                chosen_keys.add(key)

    clean_sym = re.sub(r"[^A-Za-z0-9_-]", "", gene_symbol.upper()) or "GENE"
    can_design_primers = len(seq) >= 380
    for rank, g in enumerate(selected, start=1):
        g.guide_id = f"{clean_sym}_E{exon_number}_g{rank}"
        if can_design_primers:
            g.genotyping_primers = design_tide_ice_primers(
                window_plus_seq=seq,
                window_genomic_start=p_start,
                cut_genomic_pos=g.genomic_cut_pos,
                chrom=chrom,
            )

    return selected


def design_crispri_tss_guides(
    tss_window_plus_seq: str,
    window_genomic_start: int,
    tss_pos: int,
    chrom: str,
    gene_strand: int,
    gene_symbol: str = "GENE",
    max_guides: int = 5,
    nuclease_id: str = "spcas9",
) -> list[GuideCandidate]:
    """Design dCas9/dCas12a-KRAB CRISPRi promoter guides in [TSS - 50 bp .. TSS + 300 bp]."""
    if not tss_window_plus_seq or tss_pos <= 0:
        return []
    nuc = resolve_nuclease(nuclease_id)
    seq = tss_window_plus_seq.upper()
    chr_label = chrom if chrom.startswith("chr") else f"chr{chrom.replace('chr', '')}"
    n = len(seq)

    raw_hits = _scan_nuclease_sites(seq, window_genomic_start, 0, n, nuc)
    cands: list[GuideCandidate] = []
    seen: set[tuple[str, int]] = set()

    for spacer, pam, g_strand, target_pos, g_s_win, g_e_win in raw_hits:
        if "TTTT" in spacer or (spacer, target_pos) in seen:
            continue
        tss_off = (target_pos - tss_pos) if gene_strand >= 0 else (tss_pos - target_pos)
        if tss_off < -50 or tss_off > 300:
            continue
        seen.add((spacer, target_pos))

        gc_pct, score, score_notes, fwd_oligo, rev_oligo = _score_nuclease_spacer(
            spacer=spacer,
            pam=pam,
            nuc=nuc,
            base_score=76,
        )
        notes: list[str] = []
        # Horlbeck et al. (2016) & Gilbert et al. (2014): peak dCas-KRAB repression at +25..+125 bp past TSS
        if 25 <= tss_off <= 125:
            score += 12
            notes.append(f"Peak CRISPRi window (TSS {tss_off:+d} bp)")
        elif 0 <= tss_off <= 200:
            score += 6
            notes.append(f"Strong TSS window ({tss_off:+d} bp)")
        else:
            notes.append(f"Active TSS window ({tss_off:+d} bp)")
        notes.extend(score_notes[:2])

        tx_orient = (
            "Sense"
            if ((g_strand == "+" and gene_strand >= 0) or (g_strand == "-" and gene_strand < 0))
            else "Antisense"
        )
        if nuc.nuclease_id == "cas12a":
            if tx_orient == "Sense" and 15 <= tss_off <= 150:
                score += 6  # dCas12a favors template-strand pairing (Tak et al., Nat Methods 2017)
        elif tx_orient == "Antisense" and 15 <= tss_off <= 150:
            score += 3  # dCas9 non-template strand steric roadblock bonus

        cands.append(
            GuideCandidate(
                guide_id="",
                exon_number=0,
                spacer_20nt=spacer,
                pam=pam,
                genomic_strand=g_strand,
                transcript_orientation=tx_orient,
                chrom=chr_label,
                genomic_cut_pos=target_pos,
                aa_cut=0,
                gc_percent=gc_pct,
                overlaps_domain=f"Promoter TSS {tss_off:+d} bp",
                sequence_score=max(0, min(100, score)),
                quality_notes=notes,
                cloning_fwd_oligo=fwd_oligo,
                cloning_rev_oligo=rev_oligo,
                mode="CRISPRi",
                tss_offset_bp=tss_off,
                genomic_start_23bp=g_s_win,
                genomic_end_23bp=g_e_win,
                nuclease=nuc.short_tag,
                spacer_length=nuc.spacer_length,
                pam_orientation=nuc.pam_orientation,
                cloning_backbone_hint=nuc.cloning_backbone,
            )
        )

    cands.sort(key=lambda g: (-g.sequence_score, abs((g.tss_offset_bp or 0) - 65), g.genomic_cut_pos))
    selected: list[GuideCandidate] = []
    for cand in cands:
        if len(selected) >= max_guides:
            break
        if all(abs(cand.genomic_cut_pos - prev.genomic_cut_pos) >= 12 for prev in selected):
            selected.append(cand)

    clean_sym = re.sub(r"[^A-Za-z0-9_-]", "", gene_symbol.upper()) or "GENE"
    for rank, g in enumerate(selected, start=1):
        g.guide_id = f"{clean_sym}_CRISPRi_g{rank}"
    return selected


def select_top_gene_guides(exons: list[ExonAudit], max_total: int = 6) -> list[GuideCandidate]:
    """Select a diverse multi-exon panel of top candidate sgRNAs across the safest coding exons."""
    coding_with_guides = [
        e for e in exons if e.is_coding and e.candidate_guides
    ]
    if not coding_with_guides:
        return []

    rec_exons = [e for e in coding_with_guides if e.tier == "RECOMMENDED"]
    if rec_exons:
        pool = sorted(rec_exons, key=lambda e: (-e.priority_score, e.exon_number))
    else:
        caut_exons = [e for e in coding_with_guides if e.priority_score >= 45]
        pool = sorted(
            caut_exons or coding_with_guides,
            key=lambda e: (-e.priority_score, e.exon_number),
        )

    top_exons = pool[:3]
    panel: list[GuideCandidate] = []
    round_idx = 0
    while len(panel) < max_total:
        added_in_round = False
        for ex in top_exons:
            if round_idx < len(ex.candidate_guides) and len(panel) < max_total:
                panel.append(ex.candidate_guides[round_idx])
                added_in_round = True
        if not added_in_round:
            break
        round_idx += 1

    return panel


def score_transcript_exons(
    chrom: str,
    primary_tx: dict[str, Any],
    alt_coding_txs: list[dict[str, Any]],
    domains: list[ProteinDomain],
    cds_seq: str = "",
    flanked_seqs: Optional[dict[tuple[int, int], Optional[str]]] = None,
    flank_pad: int = 25,
    gene_symbol: str = "GENE",
    nuclease_id: str = "spcas9",
) -> list[ExonAudit]:
    """Audit and score every exon of the primary transcript in 5' -> 3' biological order."""
    nuc = resolve_nuclease(nuclease_id)
    strand = int(primary_tx.get("strand", 1))
    raw_exons = list(primary_tx.get("Exon") or [])
    if strand >= 0:
        ordered_exons = sorted(raw_exons, key=lambda e: int(e["start"]))
    else:
        ordered_exons = sorted(raw_exons, key=lambda e: int(e["end"]), reverse=True)

    total_exons = len(ordered_exons)
    trans = primary_tx["Translation"]
    t_gstart = min(int(trans["start"]), int(trans["end"]))
    t_gend = max(int(trans["start"]), int(trans["end"]))

    # Pre-collect coding genomic intervals for each alternative isoform
    alt_isoform_intervals: list[tuple[str, list[tuple[int, int]]]] = []
    for alt in alt_coding_txs:
        alt_name = str(alt.get("display_name") or alt.get("id"))
        alt_tr = alt.get("Translation") or {}
        if alt_tr.get("start") is None or alt_tr.get("end") is None:
            continue
        ag_s, ag_e = min(int(alt_tr["start"]), int(alt_tr["end"])), max(int(alt_tr["start"]), int(alt_tr["end"]))
        ivals: list[tuple[int, int]] = []
        for aex in alt.get("Exon") or []:
            es, ee = min(int(aex["start"]), int(aex["end"])), max(int(aex["start"]), int(aex["end"]))
            ovs, ove = max(es, ag_s), min(ee, ag_e)
            if ovs <= ove:
                ivals.append((ovs, ove))
        alt_isoform_intervals.append((alt_name, ivals))

    # Total spliced cDNA length and final exon-exon splice junction position
    exon_lengths = [abs(int(e["end"]) - int(e["start"])) + 1 for e in ordered_exons]
    total_cdna_len = sum(exon_lengths)
    final_junction_cdna = (total_cdna_len - exon_lengths[-1]) if total_exons >= 2 else 0
    escape_zone_start_cdna = max(1, final_junction_cdna - 55 + 1) if total_exons >= 2 else 1

    # Downstream in-frame ATG codons & earliest functional domain start
    internal_mets = find_internal_inframe_mets(cds_seq, max_nt=450)
    first_domain_aa = min((d.aa_start for d in domains), default=999999)

    clean_chr = str(chrom).replace("chr", "")
    chr_label = f"chr{clean_chr}"

    audits: list[ExonAudit] = []
    cum_cdna = 0
    cum_cds = 0

    for idx, ex in enumerate(ordered_exons, start=1):
        g_start = min(int(ex["start"]), int(ex["end"]))
        g_end = max(int(ex["start"]), int(ex["end"]))
        ex_len = g_end - g_start + 1
        cdna_start = cum_cdna + 1
        cum_cdna += ex_len

        ov_start = max(g_start, t_gstart)
        ov_end = min(g_end, t_gend)
        if ov_start > ov_end:
            audits.append(
                ExonAudit(
                    exon_number=idx,
                    exon_id=str(ex.get("id") or f"Exon_{idx}"),
                    chrom=chr_label,
                    genomic_start=g_start,
                    genomic_end=g_end,
                    strand=strand,
                    exon_length_bp=ex_len,
                    is_coding=False,
                    nmd_status="NON_CODING_UTR",
                    priority_score=0,
                    tier="UTR_NONCODING",
                    score_breakdown=["Non-coding UTR exon (cutting here does not frameshift the CDS)"],
                )
            )
            continue

        coding_len = ov_end - ov_start + 1
        cds_off_start = cum_cds + 1
        cds_off_end = cum_cds + coding_len
        cum_cds += coding_len

        aa_start = ((cds_off_start - 1) // 3) + 1
        aa_end = ((cds_off_end - 1) // 3) + 1

        if strand >= 0:
            coding_cdna_start = cdna_start + (ov_start - g_start)
        else:
            coding_cdna_start = cdna_start + (g_end - ov_end)
        coding_cdna_end = coding_cdna_start + coding_len - 1

        # 1. Isoform constitutiveness across primary + well-supported coding isoforms
        skipped_in: list[str] = []
        for alt_name, ivals in alt_isoform_intervals:
            # Consider retained if >=80% of this exon's coding bp overlap an alt coding exon
            covered_bp = 0
            for a_s, a_e in ivals:
                is_s, is_e = max(ov_start, a_s), min(ov_end, a_e)
                if is_s <= is_e:
                    covered_bp += is_e - is_s + 1
            if covered_bp < 0.8 * coding_len:
                skipped_in.append(alt_name)

        total_iso_compared = 1 + len(alt_isoform_intervals)
        retained_count = total_iso_compared - len(skipped_in)
        iso_fraction = round(retained_count / total_iso_compared, 3)
        is_constitutive = len(skipped_in) == 0

        # 2. Internal exon symmetry (in-frame skipping risk when entire internal exon is coding and divisible by 3)
        is_internal_wholly_coding = (
            1 < idx < total_exons and coding_len == ex_len
        )
        is_symmetric = bool(is_internal_wholly_coding and (ex_len % 3 == 0))

        # 3. NMD escape & 5' translation re-initiation geometry
        reinit_met: Optional[int] = None
        for m_aa in internal_mets:
            if m_aa > aa_start and m_aa <= first_domain_aa + 15:
                reinit_met = m_aa
                break

        if total_exons == 1:
            nmd_escape_frac = 1.0
            nmd_status = "ESCAPE_SINGLE_EXON"
        else:
            escape_bp = max(
                0,
                coding_cdna_end - max(coding_cdna_start, escape_zone_start_cdna) + 1,
            )
            nmd_escape_frac = round(escape_bp / coding_len, 3)
            if idx == total_exons:
                nmd_status = "ESCAPE_LAST_EXON"
            elif nmd_escape_frac >= 0.5:
                nmd_status = "ESCAPE_FINAL_55NT"
            elif cds_off_start <= 150 and reinit_met is not None:
                nmd_status = "REINIT_5PRIME"
            elif ex_len > 400:
                nmd_status = "ATTENUATED_LONG_EXON"
            else:
                nmd_status = "COMPETENT"

        # 4. Conserved functional domain overlap
        hit_domains: list[str] = []
        for dom in domains:
            if max(aa_start, dom.aa_start) <= min(aa_end, dom.aa_end):
                label = f"{dom.description} ({dom.source_db}:{dom.accession})"
                if label not in hit_domains:
                    hit_domains.append(label)

        # 5. Nuclease PAM count & candidate sgRNA/crRNA extraction inside coding interval
        pam_count: Optional[int] = None
        exon_guides: list[GuideCandidate] = []
        if flanked_seqs is not None and (ov_start, ov_end) in flanked_seqs:
            fseq = flanked_seqs[(ov_start, ov_end)]
            if fseq:
                actual_left_pad = min(flank_pad, ov_start - 1)
                pam_count = count_spcas9_cut_sites(fseq, actual_left_pad, coding_len, nuclease_id=nuc.nuclease_id)
                exon_guides = design_exon_sgrnas(
                    padded_forward_seq=fseq,
                    left_pad=actual_left_pad,
                    coding_len=coding_len,
                    chrom=chr_label,
                    coding_genomic_start=ov_start,
                    coding_genomic_end=ov_end,
                    gene_strand=strand,
                    cds_offset_start=cds_off_start,
                    exon_number=idx,
                    domains=domains,
                    gene_symbol=gene_symbol,
                    is_single_exon=(total_exons == 1),
                    max_guides=5,
                    nuclease_id=nuc.nuclease_id,
                )

        # 6. Transparent additive Exon Priority Score (base = 80, clamped to 0..100)
        score = 80
        breakdown: list[str] = ["Base coding exon score: 80"]

        if hit_domains:
            score += 20
            short_dom = hit_domains[0].split(" (")[0]
            breakdown.append(
                f"+20 Encodes conserved domain '{short_dom}' (even -3 bp in-frame indels disrupt function)"
            )

        if nmd_status == "ESCAPE_SINGLE_EXON":
            if hit_domains:
                score -= 25
                breakdown.append(
                    "-25 Single-exon transcript: frameshifts escape EJC-NMD (relying on domain disruption)"
                )
            else:
                score -= 55
                breakdown.append(
                    "-55 Single-exon transcript outside annotated domains: frameshifts escape NMD"
                )
        elif nmd_status in ("ESCAPE_LAST_EXON", "ESCAPE_FINAL_55NT"):
            score -= 55
            zone_desc = (
                "last exon" if nmd_status == "ESCAPE_LAST_EXON" else "final 55-nt EJC-blind zone"
            )
            breakdown.append(
                f"-55 {int(nmd_escape_frac * 100)}% of coding span sits in {zone_desc} (NMD escape -> truncated protein)"
            )
        elif nmd_escape_frac > 0.0:
            score -= 15
            breakdown.append(
                f"-15 3' end of exon borders final 55-nt NMD escape zone ({int(nmd_escape_frac * 100)}% overlap)"
            )

        if not is_constitutive:
            penalty = 25 if iso_fraction <= 0.67 else 15
            score -= penalty
            breakdown.append(
                f"-{penalty} Spliced out in {len(skipped_in)} coding isoform(s): {', '.join(skipped_in[:3])}"
            )

        if cds_off_start <= 150 and reinit_met is not None:
            score -= 20
            breakdown.append(
                f"-20 Within first 150 nt of CDS with downstream in-frame Met{reinit_met} (translation re-initiation risk)"
            )

        if is_symmetric:
            score -= 15
            breakdown.append(
                f"-15 Symmetric internal exon ({ex_len} bp, divisible by 3 -> in-frame exon skipping risk)"
            )

        if ex_len > 400 and nmd_status not in ("ESCAPE_LAST_EXON", "ESCAPE_SINGLE_EXON"):
            score -= 15
            breakdown.append(
                f"-15 Long exon ({ex_len} bp > 400 bp -> attenuated NMD efficiency per Lindeboom 2019)"
            )

        if pam_count is not None:
            if pam_count == 0:
                score -= 45
                alt_hint = (
                    "switch Nuclease to SpG [NGN] or Cas12a [TTTV]"
                    if nuc.nuclease_id != "spg"
                    else "select an adjacent exon"
                )
                breakdown.append(
                    f"-45 Zero {nuc.short_tag} ({nuc.pam_pattern}) cut sites inside coding interval ({alt_hint})"
                )
            elif pam_count <= 2:
                score -= 10
                breakdown.append(
                    f"-10 Sparse {nuc.short_tag} ({nuc.pam_pattern}) options ({pam_count} cut site(s) in coding span)"
                )

        score = max(0, min(100, score))
        if score >= 75:
            tier = "RECOMMENDED"
        elif score >= 45:
            tier = "CAUTION"
        else:
            tier = "AVOID"

        audits.append(
            ExonAudit(
                exon_number=idx,
                exon_id=str(ex.get("id") or f"Exon_{idx}"),
                chrom=chr_label,
                genomic_start=g_start,
                genomic_end=g_end,
                strand=strand,
                exon_length_bp=ex_len,
                is_coding=True,
                coding_genomic_start=ov_start,
                coding_genomic_end=ov_end,
                coding_length_bp=coding_len,
                cds_offset_start=cds_off_start,
                cds_offset_end=cds_off_end,
                aa_start=aa_start,
                aa_end=aa_end,
                is_constitutive=is_constitutive,
                isoform_coverage_fraction=iso_fraction,
                skipped_in_isoforms=skipped_in,
                is_symmetric=is_symmetric,
                nmd_status=nmd_status,
                nmd_escape_fraction=nmd_escape_frac,
                reinitiation_met_codon=reinit_met if cds_off_start <= 150 else None,
                overlapping_domains=hit_domains,
                spcas9_pam_count=pam_count,
                priority_score=score,
                tier=tier,
                score_breakdown=breakdown,
                crispick_target_region=f"{chr_label}:{ov_start}-{ov_end}",
                candidate_guides=exon_guides,
            )
        )

    return audits
