"""CRISPR nuclease / Cas-ortholog registry, PAM matching, and cloning overhang rules.

Supports 5 wet-lab nuclease configurations using 100% Python stdlib:
- `spcas9`: Wild-type SpCas9 (`20-nt + 3'-NGG`, blunt cut 3 bp 5' of PAM)
- `spcas9_hf`: High-Fidelity SpCas9 (`eSpCas9(1.1)` / `SpCas9-HF1` / `HiFi Cas9`, `20-nt + 3'-NGG`, native `5'-G` priority)
- `spg`: PAM-relaxed SpG Cas9 (`20-nt + 3'-NGN`, Walton et al., Science 2020)
- `sacas9`: Compact AAV SaCas9 (`21-nt + 3'-NNGRRT`, Ran et al., Nature 2015; BsaI pX601)
- `cas12a`: AsCas12a / LbCas12a (`5'-TTTV + 23-nt`, Zetsche et al., Cell 2015; staggered cut & crRNA BsmBI)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class NucleaseConfig:
    """Biochemical and cloning specification for a CRISPR nuclease."""

    nuclease_id: str
    display_name: str
    short_tag: str
    pam_pattern: str
    pam_orientation: str  # "3prime" ([spacer][PAM]) or "5prime" ([PAM][spacer])
    spacer_length: int
    pam_length: int
    cloning_backbone: str

    @property
    def window_length(self) -> int:
        return self.spacer_length + self.pam_length

    @property
    def target_format_label(self) -> str:
        if self.pam_orientation == "5prime":
            return f"5'-{self.pam_pattern} + {self.spacer_length}-nt"
        return f"{self.spacer_length}-nt + 3'-{self.pam_pattern}"


NUCLEASE_REGISTRY: dict[str, NucleaseConfig] = {
    "spcas9": NucleaseConfig(
        nuclease_id="spcas9",
        display_name="SpCas9 (Standard 20-nt + NGG)",
        short_tag="SpCas9",
        pam_pattern="NGG",
        pam_orientation="3prime",
        spacer_length=20,
        pam_length=3,
        cloning_backbone="U6 BsmBI/BbsI (lentiCRISPRv2 / PX459)",
    ),
    "spcas9_hf": NucleaseConfig(
        nuclease_id="spcas9_hf",
        display_name="HiFi / eSpCas9 / HF1 (20-nt + NGG, native 5'-G)",
        short_tag="SpCas9-HF",
        pam_pattern="NGG",
        pam_orientation="3prime",
        spacer_length=20,
        pam_length=3,
        cloning_backbone="U6 BsmBI/BbsI (eSpCas9 / HF1 — native 5'-G)",
    ),
    "spg": NucleaseConfig(
        nuclease_id="spg",
        display_name="SpG Cas9 (PAM-Relaxed 20-nt + NGN)",
        short_tag="SpG",
        pam_pattern="NGN",
        pam_orientation="3prime",
        spacer_length=20,
        pam_length=3,
        cloning_backbone="U6 BsmBI/BbsI (SpG vector)",
    ),
    "sacas9": NucleaseConfig(
        nuclease_id="sacas9",
        display_name="SaCas9 (Compact AAV 21-nt + NNGRRT)",
        short_tag="SaCas9",
        pam_pattern="NNGRRT",
        pam_orientation="3prime",
        spacer_length=21,
        pam_length=6,
        cloning_backbone="U6 BsaI (AAV pX601 / pX602)",
    ),
    "cas12a": NucleaseConfig(
        nuclease_id="cas12a",
        display_name="AsCas12a / LbCas12a (5'-TTTV + 23-nt)",
        short_tag="Cas12a",
        pam_pattern="TTTV",
        pam_orientation="5prime",
        spacer_length=23,
        pam_length=4,
        cloning_backbone="U6 crRNA BsmBI (pY010 / pRDA_052)",
    ),
}

NUCLEASE_ALIASES: dict[str, str] = {
    "spcas9": "spcas9",
    "cas9": "spcas9",
    "ngg": "spcas9",
    "wt": "spcas9",
    "spcas9_hf": "spcas9_hf",
    "spcas9-hf": "spcas9_hf",
    "spcas9-hf1": "spcas9_hf",
    "hf1": "spcas9_hf",
    "espcas9": "spcas9_hf",
    "hifi": "spcas9_hf",
    "spg": "spg",
    "ngn": "spg",
    "sacas9": "sacas9",
    "sa": "sacas9",
    "aav": "sacas9",
    "nngrrt": "sacas9",
    "cas12a": "cas12a",
    "cas12": "cas12a",
    "cpf1": "cas12a",
    "ascas12a": "cas12a",
    "lbcas12a": "cas12a",
    "enascas12a": "cas12a",
    "tttv": "cas12a",
}


def resolve_nuclease(raw_nuclease: Optional[str]) -> NucleaseConfig:
    """Resolve a user-supplied nuclease identifier or alias to a NucleaseConfig."""
    if not raw_nuclease or not raw_nuclease.strip():
        return NUCLEASE_REGISTRY["spcas9"]
    key = raw_nuclease.strip().lower().replace(" ", "_")
    resolved = NUCLEASE_ALIASES.get(key)
    if not resolved:
        valid = ", ".join(NUCLEASE_REGISTRY.keys())
        raise ValueError(f"Unsupported nuclease '{raw_nuclease}'. Supported nucleases: {valid}")
    return NUCLEASE_REGISTRY[resolved]
