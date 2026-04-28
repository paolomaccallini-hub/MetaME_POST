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

## External tools

- [METAL](https://github.com/statgen/METAL) — weighted Z-score meta-analysis. The pipeline runs METAL via WSL on Windows. METAL source code (version 2011-03-25) must be downloaded from [THIS PAGE](https://csg.sph.umich.edu/abecasis/Metal/download/) and compiled. To integrate in the MetaME, adjust `path_metal_exe` in `MetaME_config.yml` for your system.
- [FUMA](https://fuma.ctglab.nl) — SNP2GENE, MAGMA gene-set/tissue/cell-type analyses. FUMA jobs are submitted manually using the config files in `FUMA/configs/`.

---

## Data

### Summary Statistics

Raw summary statistics are not included in this repository due to size and access restrictions. The script will download them automatically from the following sources:

| Database | Symbol | Cases | Controls | Trait | Regression | Ancestry | Assembly | Reference | Summary Statistics |
| :------- | :----- | -----:| --------:| :---- | :--------- | :------- | :------- | :-------- | :----------------- |
| **DecodeME** | DME_1 | 15,579 | 259,909 | CFS (CCC/IOM) | Logistic | EUR | GRCh38 | [Preprint_2025](https://www.research.ed.ac.uk/en/publications/initial-findings-from-the-decodeme-genome-wide-association-study-) | [GWAS-1](https://osf.io/rgqs3/files) |
| **Million_Veteran_Program** | MVP | 3,891 | 439,202 | PheCode_798.1 CFS | Logistic (SAIGE) | EUR | GRCh38 | [Verma_2024](https://pubmed.ncbi.nlm.nih.gov/39024449/) | [GCST90479178](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90479001-GCST90480000/GCST90479178/) |

### Gene-sets, tissue expression, and cell-type datasets

The following files must be manually downloaded from FUMA [download page](https://fuma.ctglab.nl/downloadPage) and put on the main folder.

| File | Size | Description |
|:-----|-----:|:------------|
| gtex_v8_ts_DEG.txt | 12MB | The gene-set of the differentially expressed genes in GTEx v8 across 54 tissues. This file is used for GENE2FUNC tissue specificity analysis. |
| MSigDB_20231Hs_MAGMA.txt | 24MB | MAGMA gene-set analysis file used in SNP2GENE from FUMA version 1.5.6 onwards. |
| preprocessed_scrnaseq.tar.gz | 1.9GB | This is a tar zip folder containing pre-processed scRNAseq datasets for the FUMA Cell Type module. After download, use the command "tar -xzvf" to untar. |

---

## Quick start

1. Clone the repository:
```bash
git clone https://github.com/<your-username>/MetaME_POST.git
cd MetaME
```

2. Download METAL source code from [this page](https://csg.sph.umich.edu/abecasis/Metal/download/) and compile it

3. Edit `MetaME_config.yml` to set the path to your METAL executable:
```yaml
METAL:
  path_metal_exe: "/path/to/metal"
```

4. Run the main pipeline:
```r
source("MetaME_main.R")
```

5. Submit FUMA jobs using the config files in `FUMA/configs/`. Full dataset identifiers and FUMA job parameters are documented there. FUMA output is available in this repository.

6. Once FUMA results are available, place them in `FUMA/results/` and run the post-processing script:
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

Foreground gene sets for ORA replication were constructed from the FUMA preprocessed scRNA-seq matrices using the top-decile expression proportion (TDEP) as described in ([Yao S. et al. 2024](https://www.nature.com/articles/s41467-024-55611-1)) and the
following further condition

```
E(g,c) - mean_A(g) > 1
```

where `E(g,c)` is the log2 transformed pseudocount of gene *g* in cell type *c* and `mean_A(g)` is the cross-cell-type geometric mean expression for gene *g* (on the log2 scale). This restricts each foreground to the top 10% most-expressed genes in each cell type with at least 2-fold enrichment over the cross-cell-type mean.

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
