# CRISPR Preflight

**Pre-experiment biological flight check for CRISPR knockout, CRISPRi, and target validation across Human, Mouse, Rat, Zebrafish, *Drosophila*, and *C. elegans*.**

Before you order sgRNA oligos and spend weeks cloning guides, selecting cells or F0 crispants, and running Western blots, `crispr-preflight` checks **whether standard single-gene Cas9 knockout is biologically valid** in your organism and model system, **which coding exons actually destroy protein function**, **designs ready-to-order cloning oligos + genotyping primers**, and **plans bench validation without false-clean Westerns**.

- **6 flagship research organisms supported out of the box:**
  - **Human** — *Homo sapiens* (`GRCh38`)
  - **Mouse** — *Mus musculus* (`GRCm39`)
  - **Rat** — *Rattus norvegicus* (`GRCr8`)
  - **Zebrafish** — *Danio rerio* (`GRCz11`)
  - **Fruit Fly** — *Drosophila melanogaster* (`BDGP6.46`)
  - **Roundworm** — *Caenorhabditis elegans* (`WBcel235`)
- **5 CRISPR nucleases & PAM systems supported:**
  - **SpCas9** (`20-nt + 3'-NGG`, `BsmBI`/`BbsI` U6 Golden Gate)
  - **HiFi / eSpCas9 / SpCas9-HF1** (`20-nt + 3'-NGG`, native `5'-GN19` U6 matching rule)
  - **SpG Cas9** (`20-nt + 3'-NGN`, PAM-relaxed 4x exon coverage)
  - **SaCas9** (`21-nt + 3'-NNGRRT`, single-AAV `pX601` `BsaI` packaging)
  - **AsCas12a / LbCas12a** (`5'-TTTV + 23-nt`, AT-rich genomes, staggered 5-nt cuts & `BsmBI` crRNA arrays)
- **Zero third-party dependencies:** 100% Python (`>=3.9`) standard library — installs in seconds with `pip` or `uv` with zero compiled C/Rust dependencies.
- **Zero API keys required:** Queries open public biological APIs (**EMBL-EBI Proteins & InterPro**, **NCBI RefSeq & Datasets v2**, **Ensembl**, **Alliance of Genome Resources (DIOPT)**, **Open Targets v4 / DepMap Chronos**, **cBioPortal CCLE**, **gnomAD**, **GTEx v10**, **Human Protein Atlas**, and **MyGene / MyVariant**) in parallel (~1.6–2.5s cold, <0.1s cached).
- **Zero black boxes:** Every warning, `0–100` exon priority score, and guide candidate prints its exact biological math and reference genome coordinates.

---

## Why Does This Exist? (The Wet-Lab Motivation)

Systematic proteomics and RNA-seq audits of Sanger-verified frameshift "knockout" cell lines (**Smits et al., *Nature Methods* 2019**; **Tuladhar et al., *Nature Communications* 2019**) showed that **30–50% of presumed CRISPR KO lines still produce residual or truncated functional protein**.

`crispr-preflight` runs **right before ordering oligos** to catch the **6 classic bench failure modes**:

1. **Cell-Line Endogenous Mutation & Amplification Traps (`cBioPortal` CCLE across 1,739 lines):**
   - Catches genes that are **already LoF-mutated or homozygously deleted** in your chosen cell line (e.g., *SMARCA4* `p.Q729Cfs*4` in `A549`), endogenous driver mutations (e.g., *KRAS* `p.G13D` in `HCT116`), and **high-level focal amplifications (`GISTIC +2`)** where Cas9 nuclease causes multi-DSB DNA damage toxicity (**Aguirre et al., 2016**).
2. **Isoform Skipping Trap:** Cutting an exon that is naturally spliced out in an alternative coding isoform (`MANE Plus Clinical` / `GENCODE Primary`).
3. **NMD Escape & Translation Re-initiation (`NMDetective` — Lindeboom et al., *Nature Genetics* 2019):**
   - Cutting in the **last exon** or within **55 nt upstream of the final exon–exon junction** (no downstream Exon Junction Complex to trigger Nonsense-Mediated Decay).
   - Cutting **intronless / single-exon genes** (e.g., *SOX2*, *JUN*) where frameshifts never trigger EJC-dependent NMD.
   - Cutting within the **first 150 nt of the CDS** upstream of an internal in-frame `ATG` (`Met`), enabling ribosomal re-initiation.
   - Cutting unusually **long internal exons (`>400 bp`)** with attenuated NMD efficiency.
4. **Symmetric Exon Skipping & Misleading Western Blots:** Cutting an internal exon whose length is divisible by 3 (`length % 3 == 0`). ESE-disrupting indels can skip that exon cleanly in-frame — and if your Western antibody binds that skipped exon, your blot looks 100% blank while internally deleted protein persists!
5. **Essentiality & `-3 bp` Clonal Escapers:** Expanding single-cell Cas9 KO clones for a fitness-essential gene (`DepMap Chronos <= -1.0` or `gnomAD LOEUF < 0.35`) kills true frameshifts and enriches for surviving `-3 bp` / `-6 bp` in-frame escapers — recommending **CRISPRi (`dCas9-KRAB` promoter guides)**, an **inducible degron (`dTAG`)**, or **early pooled assays** instead.
6. **Paralog Buffering (De Kegel et al., *Cell Systems* 2021):** Muted phenotype in a true KO because a high-identity human paralog (e.g., *SMARCA2* buffering *SMARCA4*) is co-expressed in the same tissue, requiring a **dual-guide KO** control.

---

## What Each Flight Check Produces

- **Perturbation Modality Verdict:** Standard Domain-Targeted Cas9 KO vs. Promoter CRISPRi / `dTAG` vs. Dual-Paralog KO vs. Switch Cell Line.
- **Interactive Exon, Pfam Domain & Antibody Epitope SVG Map:** Visual `5' -> 3'` protein architecture showing safe exons, conserved functional domains, and safe vs. risky Western blot antibody regions.
- **SpCas9 KO Spacers (`20-nt + NGG`) & Promoter CRISPRi Guides (`TSS -50..+300 bp`):** Filtered against U6 `TTTT` terminators, extreme GC, and Graf PAM-proximal inefficiency motifs, complete with `BsmBI` / `BbsI` Golden Gate overhangs (`lentiCRISPRv2` / `PX459`).
- **Automatic TIDE / Synthego ICE Genotyping PCR Primers:** Paired `430–740 bp` genomic PCR primers (SantaLucia nearest-neighbor PCR $T_m \approx 60\,^\circ\text{C}$) + recommended Sanger sequencing primer for every KO guide.
- **Common Human SNP Shield (`gnomAD` hg38 MAF >= 1%):** Flags common population variants inside the `PAM` or `Seed` region that can abolish Cas9 cleavage.
- **1-Click Exports:** IDT Bulk Order Sheet (`.csv`), standalone Benchling-ready `.html`, `.md`, and `.json`.

---

## Setup & Quickstart (From a Fresh Clone)

**Requirements:** Python `3.9+` and an internet connection (no API keys or third-party packages required).

### Option A: Standard `python3` + `venv` (works everywhere)

```bash
git clone https://github.com/MulhamNajajreh/crispr-preflight.git
cd crispr-preflight

python3 -m venv .venv
source .venv/bin/activate        # On Windows: .venv\Scripts\activate
pip install -e .
```

### Option B: Using `uv`

```bash
git clone https://github.com/MulhamNajajreh/crispr-preflight.git
cd crispr-preflight
uv sync
```

---

## Running `crispr-preflight`

### 1. Launch the Interactive Web App (Recommended)

```bash
crispr-preflight serve
```

Then open **`http://127.0.0.1:7860`** in your browser to search single genes, browse cell lines by lineage, triage 5–30 screen hits side-by-side, and export IDT oligo/primer order sheets.

### 2. Run a Single-Gene Flight Check in the Terminal

Try auditing human `SMARCA4` in lung (`A549` cells) or mouse `Smarca4` on `GRCm39`:

```bash
# Human (default: GRCh38)
crispr-preflight check SMARCA4 --tissue lung --cell-line A549 --html SMARCA4_report.html

# Mouse (GRCm39) with single-AAV SaCas9 (21-nt + NNGRRT) or C. elegans with Cas12a (5'-TTTV)
crispr-preflight check Smarca4 --species mouse --nuclease sacas9 --html Mouse_Smarca4_AAV.html
crispr-preflight check cep-1 --species worm --nuclease cas12a
```

### 3. Batch Triage a Hit List (e.g. from RNA-seq or a CRISPR Screen)

```bash
crispr-preflight batch SMARCA4,SOX2,PLK1,STAT3 --tissue lung --cell-line A549 \
  --output batch_summary.csv

# Multi-species batch triage (e.g. Mouse screen hits)
crispr-preflight batch Trp53,Smarca4,Plk1,Sox2 --species mouse --output mouse_hits.csv
```

### 4. Run with Docker (or Hugging Face Spaces)

```bash
docker build -t crispr-preflight .
docker run --rm -p 7860:7860 crispr-preflight
```

---

## Transparent Exon Priority Score (`0–100`)

Every coding exon starts at **`80`** and applies explicit, inspectable adjustments:

| Adjustment | Biological Rule | Primary Literature |
| :---: | :--- | :--- |
| **`+20`** | Overlaps a conserved **Pfam / UniProt / Smart** functional domain (even `-3 bp` in-frame indels disrupt function) | Shi et al., *Nat Biotechnol* 2015 |
| **`-55`** | Located in **Last Exon** or **Final 55-nt EJC-Blind Zone** (or intronless gene outside domains) | Lindeboom et al., *Nat Genet* 2019 |
| **`-25`** | **Non-constitutive** (spliced out in well-supported coding isoforms) | Ensembl / MANE |
| **`-20`** | **5′ Re-initiation Zone** (`<=150 nt` CDS with downstream in-frame `ATG` before first domain) | Tuladhar et al., *Nat Commun* 2019 |
| **`-15`** | **Symmetric internal exon** (`length % 3 == 0` -> in-frame skipping risk) | Smits et al., *Nat Methods* 2019 |
| **`-15`** | **Long internal exon** (`>400 bp` -> attenuated NMD efficiency) | Lindeboom et al., *Nat Genet* 2019 |
| **`-45`** | **Zero SpCas9 `NGG` cut sites** inside coding interval (`-10` if `1–2` sites) | SpCas9 PAM geometry |

Exons are classified as **`RECOMMENDED (>=75)`**, **`CAUTION (45–74)`**, or **`AVOID (<45)`**, with exact `GRCh38` coding coordinates.

---

## Running Tests

```bash
python3 -m unittest discover -s tests -v
```

---

## Data Sources & Terms of Use

`crispr-preflight` aggregates live annotations from open biological databases. When publishing or redistributing outputs, please cite and comply with each resource's terms:

- **EMBL-EBI Proteins & InterPro APIs (`Pfam` / `UniProt`):** https://www.ebi.ac.uk/about/terms-of-use
- **NCBI E-utilities & Datasets v2 (`GRCh38`, `GRCm39`, `GRCr8`, `GRCz11`, `BDGP6.46`, `WBcel235`):** https://www.ncbi.nlm.nih.gov/home/about/policies/
- **Alliance of Genome Resources (`MGI`, `RGD`, `ZFIN`, `FlyBase`, `WormBase` DIOPT Paralogs & Orthologs):** https://www.alliancegenome.org/
- **Ensembl REST API:** https://useast.ensembl.org/info/about/legal/disclaimer.html
- **Open Targets Platform GraphQL API (DepMap Chronos & Chemical Probes):** https://platform-docs.opentargets.org/licence
- **cBioPortal / CCLE (`ccle_broad_2019`):** https://www.cbioportal.org/
- **gnomAD (Genome Aggregation Database):** https://gnomad.broadinstitute.org/policies
- **GTEx Portal API v2 (`v10`):** https://gtexportal.org/home/license
- **Human Protein Atlas (HPA):** https://www.proteinatlas.org/about/licence
- **MyGene.info & MyVariant.info (`hg38` / `gnomAD`):** https://mygene.info/ & https://myvariant.info/

---

## Author & Acknowledgments

Developed by **Mulham Najajreh** ([mulhamnajajreh.com](https://mulhamnajajreh.com)).

This work was supported in part by complimentary access to Google Antigravity and Gemini 4 Argon to accelerate scientific workflows.

## License

Code is released under the [MIT License](LICENSE).
