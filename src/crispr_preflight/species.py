"""Registry of supported research model organisms, reference assemblies, and RefSeq chromosome maps."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class SpeciesConfig:
    """Configuration and reference genome metadata for a supported model organism."""

    species_id: str
    scientific_name: str
    common_name: str
    taxid: int
    ensembl_name: str
    assembly: str
    chr_prefix: str  # "chr" for vertebrates; "" or "chr" per convention
    nc_chrom_map: dict[str, str]


GRCH38_NC_MAP: dict[str, str] = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X": "NC_000023.11",
    "Y": "NC_000024.10",
    "MT": "NC_012920.1",
}

GRCM39_NC_MAP: dict[str, str] = {
    "1": "NC_000067.7",
    "2": "NC_000068.8",
    "3": "NC_000069.7",
    "4": "NC_000070.7",
    "5": "NC_000071.7",
    "6": "NC_000072.7",
    "7": "NC_000073.7",
    "8": "NC_000074.7",
    "9": "NC_000075.7",
    "10": "NC_000076.7",
    "11": "NC_000077.7",
    "12": "NC_000078.7",
    "13": "NC_000079.7",
    "14": "NC_000080.7",
    "15": "NC_000081.7",
    "16": "NC_000082.7",
    "17": "NC_000083.7",
    "18": "NC_000084.7",
    "19": "NC_000085.7",
    "X": "NC_000086.8",
    "Y": "NC_000087.8",
    "MT": "NC_005089.1",
    "M": "NC_005089.1",
}

GRCR8_NC_MAP: dict[str, str] = {
    "1": "NC_086019.1",
    "2": "NC_086020.1",
    "3": "NC_086021.1",
    "4": "NC_086022.1",
    "5": "NC_086023.1",
    "6": "NC_086024.1",
    "7": "NC_086025.1",
    "8": "NC_086026.1",
    "9": "NC_086027.1",
    "10": "NC_086028.1",
    "11": "NC_086029.1",
    "12": "NC_086030.1",
    "13": "NC_086031.1",
    "14": "NC_086032.1",
    "15": "NC_086033.1",
    "16": "NC_086034.1",
    "17": "NC_086035.1",
    "18": "NC_086036.1",
    "19": "NC_086037.1",
    "20": "NC_086038.1",
    "X": "NC_086039.1",
    "Y": "NC_086040.1",
}

GRCZ11_NC_MAP: dict[str, str] = {
    "1": "NC_007112.7",
    "2": "NC_007113.7",
    "3": "NC_007114.7",
    "4": "NC_007115.7",
    "5": "NC_007116.7",
    "6": "NC_007117.7",
    "7": "NC_007118.7",
    "8": "NC_007119.7",
    "9": "NC_007120.7",
    "10": "NC_007121.7",
    "11": "NC_007122.7",
    "12": "NC_007123.7",
    "13": "NC_007124.7",
    "14": "NC_007125.7",
    "15": "NC_007126.7",
    "16": "NC_007127.7",
    "17": "NC_007128.7",
    "18": "NC_007129.7",
    "19": "NC_007130.7",
    "20": "NC_007131.7",
    "21": "NC_007132.7",
    "22": "NC_007133.7",
    "23": "NC_007134.7",
    "24": "NC_007135.7",
    "25": "NC_007136.7",
    "MT": "NC_002333.2",
    "M": "NC_002333.2",
}

BDGP6_NC_MAP: dict[str, str] = {
    "X": "NC_004354.4",
    "2L": "NT_033779.5",
    "2R": "NT_033778.4",
    "3L": "NT_037436.4",
    "3R": "NT_033777.3",
    "4": "NC_004353.4",
    "Y": "NC_024512.1",
    "MT": "NC_024511.2",
    "M": "NC_024511.2",
}

WBCEL235_NC_MAP: dict[str, str] = {
    "I": "NC_003279.8",
    "II": "NC_003280.10",
    "III": "NC_003281.10",
    "IV": "NC_003282.8",
    "V": "NC_003283.11",
    "X": "NC_003284.9",
    "MT": "NC_001328.1",
    "M": "NC_001328.1",
}

SPECIES_REGISTRY: dict[str, SpeciesConfig] = {
    "human": SpeciesConfig(
        species_id="human",
        scientific_name="Homo sapiens",
        common_name="Human",
        taxid=9606,
        ensembl_name="homo_sapiens",
        assembly="GRCh38",
        chr_prefix="chr",
        nc_chrom_map=GRCH38_NC_MAP,
    ),
    "mouse": SpeciesConfig(
        species_id="mouse",
        scientific_name="Mus musculus",
        common_name="Mouse",
        taxid=10090,
        ensembl_name="mus_musculus",
        assembly="GRCm39",
        chr_prefix="chr",
        nc_chrom_map=GRCM39_NC_MAP,
    ),
    "rat": SpeciesConfig(
        species_id="rat",
        scientific_name="Rattus norvegicus",
        common_name="Rat",
        taxid=10116,
        ensembl_name="rattus_norvegicus",
        assembly="GRCr8",
        chr_prefix="chr",
        nc_chrom_map=GRCR8_NC_MAP,
    ),
    "zebrafish": SpeciesConfig(
        species_id="zebrafish",
        scientific_name="Danio rerio",
        common_name="Zebrafish",
        taxid=7955,
        ensembl_name="danio_rerio",
        assembly="GRCz11",
        chr_prefix="chr",
        nc_chrom_map=GRCZ11_NC_MAP,
    ),
    "fly": SpeciesConfig(
        species_id="fly",
        scientific_name="Drosophila melanogaster",
        common_name="Fruit Fly",
        taxid=7227,
        ensembl_name="drosophila_melanogaster",
        assembly="BDGP6.46",
        chr_prefix="chr",
        nc_chrom_map=BDGP6_NC_MAP,
    ),
    "worm": SpeciesConfig(
        species_id="worm",
        scientific_name="Caenorhabditis elegans",
        common_name="C. elegans",
        taxid=6239,
        ensembl_name="caenorhabditis_elegans",
        assembly="WBcel235",
        chr_prefix="chr",
        nc_chrom_map=WBCEL235_NC_MAP,
    ),
}

SPECIES_ALIASES: dict[str, str] = {
    "human": "human",
    "homo_sapiens": "human",
    "hs": "human",
    "hg38": "human",
    "grch38": "human",
    "mouse": "mouse",
    "mus_musculus": "mouse",
    "mm": "mouse",
    "mm39": "mouse",
    "grcm39": "mouse",
    "rat": "rat",
    "rattus_norvegicus": "rat",
    "rn": "rat",
    "grcr8": "rat",
    "zebrafish": "zebrafish",
    "danio_rerio": "zebrafish",
    "dr": "zebrafish",
    "zf": "zebrafish",
    "danrer11": "zebrafish",
    "grcz11": "zebrafish",
    "fly": "fly",
    "drosophila": "fly",
    "drosophila_melanogaster": "fly",
    "dm": "fly",
    "dm6": "fly",
    "bdgp6": "fly",
    "worm": "worm",
    "celegans": "worm",
    "c_elegans": "worm",
    "caenorhabditis_elegans": "worm",
    "ce": "worm",
    "ce11": "worm",
    "wbcel235": "worm",
}


def resolve_species(raw: Optional[str]) -> SpeciesConfig:
    """Normalize a user-provided species string to a supported SpeciesConfig."""
    if not raw or not raw.strip():
        return SPECIES_REGISTRY["human"]
    key = raw.strip().lower().replace(" ", "_").replace("-", "_").replace(".", "")
    mapped = SPECIES_ALIASES.get(key)
    if mapped and mapped in SPECIES_REGISTRY:
        return SPECIES_REGISTRY[mapped]
    valid = ", ".join(SPECIES_REGISTRY.keys())
    raise ValueError(f"Unsupported species '{raw}'. Supported organisms: {valid}.")
