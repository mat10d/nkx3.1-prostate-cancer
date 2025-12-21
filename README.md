# NKX3.1 Localization to Mitochondria Suppresses Prostate Cancer Initiation

This repository contains the statistical analysis code for the published research article:

**Citation:**
> Kedika VL, DiBernardo M, Sirianni J, et al. NKX3.1 Localization to Mitochondria Suppresses Prostate Cancer Initiation. *Cancer Discovery*. 2021;11(9):2316-2333. doi:10.1158/2159-8290.CD-20-1765

**Publication Link:** https://aacrjournals.org/cancerdiscovery/article-abstract/11/9/2316/666307

**Principal Investigator:** Cory Abate-Shen, PhD
**Analysis Conducted By:** Matteo Di Bernardo

---

## Overview

This repository provides reproducible R code for the clinical and genomic analyses presented in the manuscript. The analyses demonstrate that NKX3.1 localization to mitochondria (via HSPA9) suppresses prostate cancer initiation and progression.

### Key Findings Reproduced

1. **Survival Analysis**: Association between NKX3.1/HSPA9 expression and biochemical recurrence (BCR) and overall survival in multiple cohorts
2. **Co-expression Analysis**: Correlation between NKX3.1 and HSPA9 expression in prostate cancer datasets

---

## Repository Structure

```
git_repo/
├── README.md                 # This file
├── data/                     # Input data files
│   ├── HICCC.xlsx            # Included
│   ├── RP67.xlsx             # Included
│   ├── TCGA-processed-data.rds        # Download required (see below)
│   ├── TCGA-PRAD-GeneExpression-Counts.rda  # Download required (see below)
│   └── taylor-2010-processed-data.rds  # Download required (see below)
├── scripts/                  # Analysis scripts
│   ├── 01_TCGA_survival_analysis.R
│   ├── 02_TMA_HICCC_BCR_analysis.R
│   ├── 03_TMA_HICCC_overall_analysis.R
│   ├── 04_TMA_RP67_BCR_analysis.R
│   ├── 05_TCGA_coexpression_analysis.R
│   └── 06_Taylor_coexpression_analysis.R
└── results/                  # Output directory (created by scripts)
    ├── TCGA/
    ├── TMA-HICCC/
    ├── TMA-RP67/
    └── Co-expression/
```

---

## Data Download Instructions

Due to file size limitations, large data files are not included in this repository. Please download or generate the following files and place them in the `data/` directory:

### Required Data Files

| File | Size | Source |
|------|------|--------|
| `TCGA-processed-data.rds` | 82 MB | Processed TCGA PRAD clinical and expression data |
| `TCGA-PRAD-GeneExpression-Counts.rda` | 39 MB | TCGA PRAD RNA-seq counts |
| `taylor-2010-processed-data.rds` | 423 MB | Processed Taylor et al. 2010 dataset |

### Data Sources

1. **TCGA PRAD Data**: Download from [GDC Data Portal](https://portal.gdc.cancer.gov/) or [cBioPortal](https://www.cbioportal.org/study/summary?id=prad_tcga)
   - Project: TCGA-PRAD
   - Data type: Gene Expression Quantification (HTSeq Counts)
   - Clinical data: BCR status and follow-up

2. **Taylor 2010 Dataset**: Download from [GEO (GSE21035)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21035)
   - Platform: Agilent Human Genome Microarray
   - Reference: Taylor BS et al. Cancer Cell 18(1):11-22 (2010)

Contact the repository maintainer if you need assistance obtaining the processed data files.

---

## Analysis Scripts

### Survival Analysis

#### 1. TCGA Survival Analysis (`01_TCGA_survival_analysis.R`)
- **Dataset**: The Cancer Genome Atlas (TCGA) Prostate Adenocarcinoma (PRAD)
- **Analyses**:
  - Combined HSPA9 and NKX3.1 stratification (median and quartile splits)
  - HSPA9 extreme quartiles analysis
  - Biochemical recurrence-free survival
  - Cox proportional hazards models
  - Kaplan-Meier plots
- **Output**:
  - Survival plots (PDF)
  - Statistical test results (TXT)
  - Summary tables (CSV)

#### 2. TMA HICCC BCR Analysis (`02_TMA_HICCC_BCR_analysis.R`)
- **Dataset**: Herbert Irving Comprehensive Cancer Center (HICCC) Tissue Microarray
- **Analyses**:
  - Cytoplasmic, nuclear, and total NKX3.1 expression
  - Biochemical recurrence-free survival
  - Cox regression and log-rank tests
- **Output**: Kaplan-Meier plots, statistics, summary tables

#### 3. TMA HICCC Overall Survival Analysis (`03_TMA_HICCC_overall_analysis.R`)
- **Dataset**: HICCC Tissue Microarray
- **Analyses**:
  - Overall survival by NKX3.1 localization
  - Combined nuclear and cytoplasmic analysis
  - Multivariate analysis adjusting for clinical covariates:
    - Gleason score
    - T-stage
    - Age at radical prostatectomy
    - Seminal vesicle invasion (SVI)
    - PSA levels
  - Univariate and multivariate Cox models
- **Output**: Survival plots, multivariate statistics, univariate results tables

#### 4. TMA RP67 BCR Analysis (`04_TMA_RP67_BCR_analysis.R`)
- **Dataset**: RP67 Tissue Microarray
- **Analyses**:
  - Biochemical recurrence-free survival by NKX3.1 localization
  - Combined nuclear and cytoplasmic analysis
- **Output**: Kaplan-Meier plots, statistics, summary tables

### Co-expression Analysis

#### 5. TCGA Co-expression Analysis (`05_TCGA_coexpression_analysis.R`)
- **Dataset**: TCGA PRAD
- **Analyses**:
  - Pearson correlation between NKX3.1 and HSPA9
  - Analysis of both raw counts and variance-stabilized normalized expression
- **Output**: Scatter plots with regression lines, correlation statistics

#### 6. Taylor Co-expression Analysis (`06_Taylor_coexpression_analysis.R`)
- **Dataset**: Taylor et al. 2010 (GSE21035)
- **Reference**: Taylor BS et al. Cancer Cell 18(1):11-22 (2010)
- **Analyses**:
  - Correlation analysis for three HSPA9B probe sets
  - Pearson correlation with NKX3.1
- **Output**: Scatter plots for each probe, correlation statistics

---

## Requirements

### R Version
- R version 4.0.0 or higher recommended

### R Packages

Install required packages by running:

```r
# CRAN packages
install.packages(c(
  "survival",
  "survminer",
  "dplyr",
  "tidyr",
  "readxl",
  "ggplot2",
  "cowplot",
  "grid",
  "ggpubr",
  "data.table"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "Biobase"
))
```

### Complete Dependency List

**Data Manipulation:**
- dplyr
- tidyr
- data.table
- readxl

**Survival Analysis:**
- survival
- survminer

**Genomics:**
- DESeq2
- Biobase

**Visualization:**
- ggplot2
- ggpubr
- cowplot
- grid

---

## Running the Analyses

### Setting Up

1. Clone or download this repository
2. Ensure all required R packages are installed
3. Set your working directory to the `git_repo` folder

### Running Individual Scripts

Navigate to the `scripts/` directory and run scripts in order:

```r
# Set working directory to scripts folder
setwd("path/to/git_repo/scripts")

# Run survival analyses
source("01_TCGA_survival_analysis.R")
source("02_TMA_HICCC_BCR_analysis.R")
source("03_TMA_HICCC_overall_analysis.R")
source("04_TMA_RP67_BCR_analysis.R")

# Run co-expression analyses
source("05_TCGA_coexpression_analysis.R")
source("06_Taylor_coexpression_analysis.R")
```

### Running All Analyses

To run all analyses in sequence:

```r
setwd("path/to/git_repo/scripts")

scripts <- c(
  "01_TCGA_survival_analysis.R",
  "02_TMA_HICCC_BCR_analysis.R",
  "03_TMA_HICCC_overall_analysis.R",
  "04_TMA_RP67_BCR_analysis.R",
  "05_TCGA_coexpression_analysis.R",
  "06_Taylor_coexpression_analysis.R"
)

for (script in scripts) {
  cat("\n========================================\n")
  cat("Running:", script, "\n")
  cat("========================================\n\n")
  source(script)
}
```

---

## Output Files

### Results Directory Structure

After running all scripts, the `results/` directory will contain:

```
results/
├── TCGA/
│   ├── HSPA9_NKX3.1_bcr_combined.pdf
│   ├── HSPA9_NKX3.1_bcr_combined_statistics.txt
│   ├── HSPA9_NKX3.1_sample_summary.csv
│   ├── HSPA9_extreme_bcr.pdf
│   ├── HSPA9_extreme_bcr_statistics.txt
│   └── HSPA9_extreme_sample_summary.csv
├── TMA-HICCC/
│   ├── BCR/
│   │   ├── cyto_BCR_KM_plot.pdf
│   │   ├── cyto_BCR_statistics.txt
│   │   ├── cyto_BCR_summary.csv
│   │   ├── nucl_BCR_KM_plot.pdf
│   │   ├── nucl_BCR_statistics.txt
│   │   ├── nucl_BCR_summary.csv
│   │   ├── total_BCR_KM_plot.pdf
│   │   ├── total_BCR_statistics.txt
│   │   └── total_BCR_summary.csv
│   └── Overall/
│       ├── cyto_overall_KM_plot.pdf
│       ├── cyto_overall_statistics.txt
│       ├── cyto_overall_summary.csv
│       ├── cyto_overall_multivariate_statistics.txt
│       ├── cyto_overall_univariate_results.csv
│       ├── nucl_overall_KM_plot.pdf
│       ├── nucl_overall_statistics.txt
│       ├── nucl_overall_summary.csv
│       ├── nucl_overall_multivariate_statistics.txt
│       ├── nucl_overall_univariate_results.csv
│       ├── total_overall_KM_plot.pdf
│       ├── total_overall_statistics.txt
│       ├── total_overall_summary.csv
│       ├── total_overall_multivariate_statistics.txt
│       ├── total_overall_univariate_results.csv
│       ├── nuc_cyto_overall_combined_KM_plot.pdf
│       └── nuc_cyto_overall_statistics.txt
├── TMA-RP67/
│   └── BCR/
│       ├── cyto_BCR_KM_plot.pdf
│       ├── cyto_BCR_statistics.txt
│       ├── cyto_BCR_summary.csv
│       ├── nucl_BCR_KM_plot.pdf
│       ├── nucl_BCR_statistics.txt
│       ├── nucl_BCR_summary.csv
│       ├── total_BCR_KM_plot.pdf
│       ├── total_BCR_statistics.txt
│       ├── total_BCR_summary.csv
│       ├── nuc_cyto_bcr_combined_KM_plot.pdf
│       ├── nuc_cyto_bcr_statistics.txt
│       └── nuc_cyto_bcr_summary.csv
└── Co-expression/
    ├── TCGA_correlation_normalized.pdf
    ├── TCGA_correlation_normalized_statistics.txt
    ├── TCGA_correlation_raw.pdf
    ├── TCGA_correlation_raw_statistics.txt
    ├── TCGA_correlation_summary.csv
    ├── Taylor_correlation_A_14_P128569.pdf
    ├── Taylor_correlation_A_14_P128569_statistics.txt
    ├── Taylor_correlation_A_14_P201971.pdf
    ├── Taylor_correlation_A_14_P201971_statistics.txt
    ├── Taylor_correlation_A_16_P37383218.pdf
    ├── Taylor_correlation_A_16_P37383218_statistics.txt
    └── Taylor_correlation_summary.csv
```

### File Types

- **PDF**: Kaplan-Meier survival plots and correlation scatter plots
- **TXT**: Detailed statistical output (Cox regression, log-rank tests, correlation tests)
- **CSV**: Summary tables with sample counts, event rates, and statistical results

---

## Data Description

### TCGA PRAD Dataset
- **Source**: The Cancer Genome Atlas
- **Data Type**: RNA-seq gene expression counts and clinical metadata
- **Samples**: Primary prostate tumors
- **Genes Analyzed**:
  - NKX3.1 (ENSG00000113013)
  - HSPA9 (ENSG00000167034)

### Tissue Microarray (TMA) Datasets

#### HICCC TMA
- **Source**: Herbert Irving Comprehensive Cancer Center
- **Data Type**: Immunohistochemistry scoring for NKX3.1 localization
- **Variables**:
  - Nuclear NKX3.1 expression (high/low)
  - Cytoplasmic NKX3.1 expression (high/low)
  - Total NKX3.1 expression (high/low)
  - Clinical outcomes: BCR-free survival, overall survival
  - Covariates: Gleason score, T-stage, age, SVI, PSA

#### RP67 TMA
- **Data Type**: Immunohistochemistry scoring for NKX3.1 localization
- **Variables**: Similar to HICCC TMA
- **Outcome**: BCR-free survival

### Taylor 2010 Dataset
- **Source**: GEO (GSE21035)
- **Platform**: Agilent microarray
- **Samples**: Prostate tumor and normal tissue
- **Probe Sets**:
  - NKX3-1: 1 probe
  - HSPA9B: 3 probes (A_14_P128569, A_14_P201971, A_16_P37383218)

---

## Statistical Methods

### Survival Analysis
- **Kaplan-Meier Method**: Non-parametric estimation of survival curves
- **Log-Rank Test**: Comparison of survival curves between groups
- **Cox Proportional Hazards Model**: Assessment of hazard ratios with confidence intervals
- **Multivariate Analysis**: Adjustment for clinical covariates

### Gene Expression Analysis
- **Variance Stabilizing Transformation (VST)**: Normalization of RNA-seq count data using DESeq2
- **Pearson Correlation**: Linear correlation between continuous gene expression values

### Expression Stratification
- **Median Split**: Samples divided into high/low expression groups
- **Quartile Split**: Top 25% (high) vs. bottom 25% (low), middle 50% excluded

---

## Reproducibility Notes

1. **Random Seeds**: Not applicable for these deterministic analyses
2. **Platform**: Analyses were performed on R version 4.x
3. **Data Processing**: All preprocessing steps are documented in script comments
4. **Path Handling**: Scripts use relative paths (`../data/`, `../results/`) for portability

---

## Contact and Support

For questions about the analysis code:
- Repository Maintainer: Matteo Di Bernardo

For questions about the research:
- Principal Investigator: Cory Abate-Shen, PhD
- Department of Molecular Pharmacology & Therapeutics
- Columbia University Irving Medical Center

---

## License

This code is provided for research and educational purposes. Please cite the original publication when using this code or data.

---

## Acknowledgments

This work was supported by research described in the original publication. We thank all collaborators and funding agencies that made this research possible.

---

**Last Updated:** December 2025
