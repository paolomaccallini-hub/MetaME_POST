# MetaME

**MetaME_POST** is the analysis pipeline accompanying the paper:

> xxx

This repository is self-contained and includes the full MetaME pipeline
(`MetaME_main.R`, `MetaME_func.R`) configured for the DME-1 and MVP cohorts,
as well as the ME/CFS-specific downstream analyses (`MetaME_POST.R`) and
all FUMA configuration files used in this study. Note that `MetaME_main.R`
and `MetaME_func.R` are general-purpose tools that support additional cohorts
(UK Biobank Neale Lab, UK Biobank EIB, Long Covid) beyond those used in this
paper; only the DME-1 and MVP cohorts were analysed here, as specified in
`Meta_analyses.csv`.

---

## Overview

The pipeline consists of two main scripts:

| Script | Description |
|---|---|
| `MetaME_main.R` | Downloads, filters, munges, and lifts over summary statistics; runs METAL meta-analysis; plots Z-scores for top loci |
| `MetaME_func.R` | Helper functions called by `MetaME_main.R` |
| `MetaME_POST.R` | Post-FUMA analyses: ORA replication against Zhang et al. (2025), cell-type gene set construction, supplementary table generation |
| `MetaME_config.yml` | Configuration file (file paths, QC thresholds) |
| `Meta_analyses.csv` | Specifies which cohorts enter each meta-analysis |

---

## Cohorts

| Cohort | Cases | Controls | Ancestry | Assembly |
|---|---|---|---|---|
| DME-1 (DecodeME) | 15,579 | 259,909 | EUR | GRCh38 |
| MVP | 3,891 | 439,202 | EUR | GRCh38 |

**Meta-analysis**: DME-1 + MVP, N_eff = 74,219 and 8,859,361 variants, METAL SAMPLESIZE scheme.

---

## Dependencies

### R packages

```r
# CRAN
install.packages(c("httr", "R.utils", "data.table", "yaml",
                   "ggplot2", "patchwork", "readxl", "stringr",
                   "openxlsx", "archive"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "MungeSumstats",
  "BSgenome.Hsapiens.NCBI.GRCh38",
  "SNPlocs.Hsapiens.dbSNP155.GRCh38",
  "SNPlocs.Hsapiens.dbSNP155.GRCh37",
  "BSgenome.Hsapiens.1000genomes.hs37d5",
  "TissueEnrich",
  "org.Hs.eg.db"
))
```

### External tools

- [METAL](https://github.com/statgen/METAL) — weighted Z-score meta-analysis. The pipeline runs METAL via WSL on Windows; adjust `path_metal_exe` in `MetaME_config.yml` for your system.
- [FUMA](https://fuma.ctglab.nl) — SNP2GENE, MAGMA gene-set/tissue/cell-type analyses. FUMA jobs are submitted manually using the config files in `FUMA/configs/`.

---

## Data

Raw summary statistics are not included in this repository due to size and access restrictions. Download them from the following sources before running the pipeline:

| Database | Symbol | Cases | Controls | Trait | Regression | Ancestry | Assembly | Reference | Summary Statistics |
| :------- | :----- | -----:| --------:| :---- | :--------- | :------- | :------- | :-------- | :----------------- |
| **DecodeME** | DME_1 | 15,579 | 259,909 | CFS (CCC/IOM) | Logistic | EUR | GRCh38 | [Preprint_2025](https://www.research.ed.ac.uk/en/publications/initial-findings-from-the-decodeme-genome-wide-association-study-) | [GWAS-1](https://osf.io/rgqs3/files) |
| **Million_Veteran_Program** | MVP | 3,891 | 439,202 | PheCode_798.1 CFS | Logistic (SAIGE) | EUR | GRCh38 | [Verma_2024](https://pubmed.ncbi.nlm.nih.gov/39024449/) | [GCST90479178](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479178/) |

The script will attempt to download DecodeME data automatically. MVP data requires dbGaP access and must be placed in `Data/MVP/` manually.

The FUMA preprocessed scRNA-seq archive (`preprocessed_scrnaseq.tar.gz`) must be
downloaded manually from the FUMA website and placed in `Data/`. To obtain it:

1. Register for a free account at [https://fuma.ctglab.nl](https://fuma.ctglab.nl)
2. Log in and navigate to **Cell type** → **About**
3. Download the file **preprocessed_scrnaseq.tar.gz** (approximately 4 GB)
4. Place it at `Data/preprocessed_scrnaseq.tar.gz`

This file contains the preprocessed single-nucleus RNA-seq expression matrices
from the Siletti et al. (2023), Seeker et al. (2023), and DropViz (Saunders et al.
2018) atlases, preprocessed by the FUMA team for use with the MAGMA cell-type
analysis pipeline. The `MetaME_POST.R` script reads this archive directly without
unpacking it, using the `archive` R package.

---

## Quick start

1. Clone the repository:
```bash
git clone https://github.com/<your-username>/MetaME_POST.git
cd MetaME
```

2. Edit `MetaME_config.yml` to set the path to your METAL executable:
```yaml
METAL:
  path_metal_exe: "/path/to/metal"
```

3. Run the main pipeline:
```r
source("MetaME_main.R")
```

4. Submit FUMA jobs using the config files in `FUMA/configs/`. Full dataset identifiers and FUMA job parameters are documented there. FUMA output is available in this repository.

5. Once FUMA results are available, place them in `FUMA/results/` and run the post-processing script:
```r
source("MetaME_POST.R")
```

---

## FUMA configuration

Nine FUMA jobs were submitted, corresponding to three GWAS inputs × three analysis types:

| Job | Input | Analysis |
|---|---|---|
| DME_1 | DME-1 GRCh37 | SNP2GENE + MAGMA |
| MVP | MVP GRCh37 | SNP2GENE + MAGMA |
| DME_1_MVP | Meta-analysis GRCh37 | SNP2GENE + MAGMA |
| DME_1_DropViz_L2 | DME-1 MAGMA genes | Cell-type (DropViz) |
| MVP_DropViz_L2 | MVP MAGMA genes | Cell-type (DropViz) |
| DME_1_MVP_DropViz_L2 | Meta MAGMA genes | Cell-type (DropViz) |
| DME_1_Siletti_Seeker_L2 | DME-1 MAGMA genes | Cell-type (Siletti/Seeker) |
| MVP_Siletti_Seeker_L2 | MVP MAGMA genes | Cell-type (Siletti/Seeker) |
| DME_1_MVP_Siletti_Seeker_L2 | Meta MAGMA genes | Cell-type (Siletti/Seeker) |

All jobs used GRCh37, EUR reference panel (1000 Genomes Phase 3), MHC region excluded (`exMHC=1`), and `magma_window=0`. The Seeker white-matter datasets were restricted to young donors (age 30–45 years).

---

## Cell-type gene set construction

Foreground gene sets for ORA replication were constructed from the FUMA preprocessed scRNA-seq matrices. A gene *g* was included in the foreground for cell type *c* if:

```
E(g,c) > Q90(c)   AND   E(g,c) - mean_A(g) > 1
```

where `Q90(c)` is the 90th percentile of expression within cell type *c* and `mean_A(g)` is the cross-cell-type geometric mean expression for gene *g* (both on the log2 scale as provided by FUMA). This restricts each foreground to the top 10% most expressed genes in each cell type with at least 2-fold enrichment over the cross-cell-type mean.

---

## Replication

Zhang et al. (2025) candidate genes (115 genes, q < 0.02) were used as an independent replication set. Over-representation analysis (ORA) was performed using a one-sided hypergeometric test against the STRING network background (17,759 genes). Bonferroni correction for replication was applied using k equal to the number of discoveries in the meta-analysis that survived the primary correction threshold.

> Zhang S. et al. *Dissecting the genetic complexity of myalgic encephalomyelitis/chronic fatigue syndrome via deep learning-powered genome analysis.* medRxiv 2025. doi:10.1101/2025.04.15.25325899

---

## Output

The pipeline produces:

- `Munged/` — QC-filtered and munged summary statistics (GRCh38 and GRCh37)
- `Output/` — METAL meta-analysis output, munged meta-analysis, Z-score plots
- `Replication/` — ORA replication results for gene sets, tissues, and cell types
- `Supplementary_Tables.xlsx` — all supplementary tables for the paper

---

## License

MIT License. See `LICENSE` for details.

---

## Citation

If you use this pipeline, please cite:

> Maccallini P. xxx

---

## Acknowledgements

This work is dedicated to the memory of Pierluigi Maccallini, whose intellectual guidance made this project possible.
