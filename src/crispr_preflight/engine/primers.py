"""Genomic PCR + Sanger sequencing primer designer for TIDE / Synthego ICE indel verification.

Literature basis:
- SantaLucia, J. Jr. (1998), PNAS 95:1460-1465: Unified nearest-neighbor DNA duplex thermodynamics.
- Brinkman et al. (2014), Nucleic Acids Research 42:e168 [TIDE] & Synthego ICE guidelines:
  450-700 bp total amplicon with sequencing primer placed ~185-295 bp from the Cas9 cut site
  and >=180 bp downstream trace window for accurate indel decomposition.
"""

from __future__ import annotations

import math
from typing import Optional

from crispr_preflight.models import GenotypingPrimerPair

_DNA_COMP = str.maketrans("ACGT", "TGCA")

# SantaLucia (1998) unified nearest-neighbor parameters: (dH in kcal/mol, dS in cal/(mol*K))
_NN_PARAMS: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2),
    "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7),
    "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4),
    "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0),
    "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2),
    "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
    "CC": (-8.0, -19.9),
}


def _revcomp(seq: str) -> str:
    return seq.upper().translate(_DNA_COMP)[::-1]


def calc_primer_tm(
    seq: str,
    monovalent_mm: float = 180.0,
    oligo_nm: float = 500.0,
) -> float:
    """Compute SantaLucia (1998) nearest-neighbor melting temperature (in °C) under standard PCR salt conditions.

    Default `monovalent_mm=180.0` reflects standard PCR buffer (50 mM KCl + 1.5 mM MgCl2 + 0.2 mM dNTPs)
    via the von Ahsen / Primer3 effective monovalent cation equivalence ([Mon+] + 120*sqrt([Mg2+]-[dNTP])).
    """
    s = seq.strip().upper()
    if len(s) < 10 or any(ch not in "ACGT" for ch in s):
        return 0.0

    # Terminal initiation parameters (SantaLucia 1998)
    dh = 0.0
    ds = 0.0
    for end_base in (s[0], s[-1]):
        if end_base in ("G", "C"):
            dh += 0.1
            ds += -2.8
        else:
            dh += 2.3
            ds += 4.1

    for i in range(len(s) - 1):
        step_dh, step_ds = _NN_PARAMS[s[i : i + 2]]
        dh += step_dh
        ds += step_ds

    # Monovalent salt correction on entropy: ds_salt = ds + 0.368 * (N - 1) * ln([Na+])
    na_molar = max(1e-4, monovalent_mm / 1000.0)
    ds_corr = ds + 0.368 * (len(s) - 1) * math.log(na_molar)

    # Non-self-complementary duplex concentration CT / 4
    ct_molar = max(1e-9, oligo_nm * 1e-9)
    r_gas = 1.987  # cal/(mol*K)
    denom = ds_corr + r_gas * math.log(ct_molar / 4.0)
    if denom >= 0:
        return 0.0
    tm_k = (dh * 1000.0) / denom
    return round(tm_k - 273.15, 1)


def _score_single_primer(seq: str) -> Optional[tuple[float, float]]:
    """Return (quality_penalty, tm_c) for a candidate 19-22 nt primer, or None if disqualified."""
    n = len(seq)
    if any(ch not in "ACGT" for ch in seq):
        return None
    if any(hp in seq for hp in ("AAAA", "TTTT", "GGGG", "CCCC")):
        return None

    gc_count = seq.count("G") + seq.count("C")
    gc_pct = (gc_count * 100.0) / n
    if gc_pct < 38.0 or gc_pct > 64.0:
        return None

    # Require 1-2 G/C bases in last 2 nucleotides (3' GC clamp), avoiding 3' GGG/CCC runs
    last2_gc = (1 if seq[-1] in "GC" else 0) + (1 if seq[-2] in "GC" else 0)
    if last2_gc == 0:
        return None
    if seq[-3:] in ("GGG", "CCC"):
        return None

    # Avoid 3' self-complementary 4-mer palindrome (primer-dimer seed)
    if seq[-4:] == _revcomp(seq[-4:]):
        return None

    tm = calc_primer_tm(seq)
    if tm < 57.2 or tm > 63.2:
        return None

    # Lower penalty is better (ideal Tm = 60.0 °C, GC = 50%, 3' terminal G/C)
    penalty = abs(tm - 60.0) * 2.2 + abs(gc_pct - 50.0) * 0.15
    if seq[-1] not in "GC":
        penalty += 1.2
    return penalty, tm


def design_tide_ice_primers(
    window_plus_seq: str,
    window_genomic_start: int,
    cut_genomic_pos: int,
    chrom: str,
) -> Optional[GenotypingPrimerPair]:
    """Design a genomic PCR + Sanger sequencing primer pair flanking a Cas9 cut site."""
    seq = window_plus_seq.upper()
    win_len = len(seq)
    cut_idx = cut_genomic_pos - window_genomic_start
    if cut_idx < 175 or (win_len - cut_idx) < 175:
        return None

    # Collect top upstream Forward primer candidates (5' start in [cut_idx - 380 .. cut_idx - 170])
    fwd_cands: list[tuple[float, int, int, str, float]] = []
    f_min = max(0, cut_idx - 380)
    f_max = max(0, cut_idx - 170)
    for pos in range(f_min, f_max, 2):
        for plen in (20, 21, 19, 22):
            if pos + plen >= cut_idx - 140:
                continue
            pseq = seq[pos : pos + plen]
            scored = _score_single_primer(pseq)
            if scored is not None:
                pen, tm = scored
                fwd_cands.append((pen, pos, pos + plen - 1, pseq, tm))
                break

    if not fwd_cands:
        return None
    fwd_cands.sort(key=lambda x: x[0])
    fwd_top = fwd_cands[:12]

    # Collect top downstream Reverse primer candidates (3' genomic end in [cut_idx + 170 .. cut_idx + 380])
    rev_cands: list[tuple[float, int, int, str, float]] = []
    r_min = min(win_len - 1, cut_idx + 170)
    r_max = min(win_len - 1, cut_idx + 380)
    for end_pos in range(r_min, r_max, 2):
        for plen in (20, 21, 19, 22):
            start_pos = end_pos - plen + 1
            if start_pos <= cut_idx + 140 or start_pos < 0:
                continue
            rseq = _revcomp(seq[start_pos : end_pos + 1])
            scored = _score_single_primer(rseq)
            if scored is not None:
                pen, tm = scored
                rev_cands.append((pen, start_pos, end_pos, rseq, tm))
                break

    if not rev_cands:
        return None
    rev_cands.sort(key=lambda x: x[0])
    rev_top = rev_cands[:12]

    best_pair: Optional[tuple[float, GenotypingPrimerPair]] = None
    chr_clean = chrom if chrom.startswith("chr") else f"chr{chrom}"

    for f_pen, f_s, _f_e, f_seq, f_tm in fwd_top:
        for r_pen, _r_s, r_e, r_seq, r_tm in rev_top:
            amp_bp = r_e - f_s + 1
            if amp_bp < 430 or amp_bp > 740:
                continue
            dt_tm = abs(f_tm - r_tm)
            if dt_tm > 2.8:
                continue
            # Avoid 3' cross-dimer between Fwd and Rev
            if f_seq[-4:] == _revcomp(r_seq[-4:]):
                continue

            off_fwd = cut_idx - f_s
            off_rev = r_e - cut_idx

            # TIDE / ICE sweet spot: sequencing primer 185-290 bp from cut site, >=180 bp trace on opposite side
            fwd_seq_dist_err = abs(off_fwd - 230)
            rev_seq_dist_err = abs(off_rev - 230)
            if fwd_seq_dist_err <= rev_seq_dist_err:
                seq_choice = "Forward"
                dist_pen = max(0.0, fwd_seq_dist_err - 45.0) * 0.08
            else:
                seq_choice = "Reverse"
                dist_pen = max(0.0, rev_seq_dist_err - 45.0) * 0.08

            amp_pen = abs(amp_bp - 560) * 0.015
            pair_score = f_pen + r_pen + dt_tm * 2.0 + dist_pen + amp_pen

            g_amp_s = window_genomic_start + f_s
            g_amp_e = window_genomic_start + r_e
            cand_obj = GenotypingPrimerPair(
                fwd_seq=f_seq,
                rev_seq=r_seq,
                fwd_tm=f_tm,
                rev_tm=r_tm,
                amplicon_bp=amp_bp,
                cut_offset_from_fwd=off_fwd,
                cut_offset_from_rev=off_rev,
                recommended_seq_primer=seq_choice,
                genomic_region=f"{chr_clean}:{g_amp_s}-{g_amp_e}",
            )
            if best_pair is None or pair_score < best_pair[0]:
                best_pair = (pair_score, cand_obj)

    return best_pair[1] if best_pair is not None else None
