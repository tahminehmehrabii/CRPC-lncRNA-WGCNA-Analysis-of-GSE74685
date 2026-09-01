# CRPC lncRNA WGCNA Analysis of GSE74685

This repository provides an R-based workflow for differential-expression analysis and Weighted Gene Co-expression Network Analysis (WGCNA) of the `GSE74685` dataset.

The workflow was developed in connection with a study investigating dysregulated long non-coding RNAs, metastatic progression, and epithelial–mesenchymal transition in castration-resistant prostate cancer.

## Reference Article

- Authors: Mehrabi T, Heidarzadehpilehrood R, Mobasheri M, Sobati T, Heshmati M, and Pirhoushiaran M
- Title: Dysregulated key long non-coding RNAs TP53TG1, RFPL1S, DLEU1, and HCG4 associated with epithelial-mesenchymal transition (EMT) in castration-resistant prostate cancer
- Journal: Advances in Cancer Biology - Metastasis
- Year: 2025
- Volume: 13
- Article number: 100132
- DOI: [10.1016/j.adcanc.2025.100132](https://doi.org/10.1016/j.adcanc.2025.100132)

## Scientific Background

Castration-resistant prostate cancer is an advanced form of prostate cancer associated with disease progression, treatment resistance, and a high risk of metastasis.

Long non-coding RNAs may contribute to cancer progression by influencing:

- gene-expression regulation;
- epithelial–mesenchymal transition;
- cellular plasticity;
- tumor-cell migration;
- invasion;
- metastatic progression;
- treatment resistance.

The reference article analyzed the `GSE74685` dataset and reported four dysregulated lncRNAs:

- `TP53TG1`
- `RFPL1S`
- `DLEU1`
- `HCG4`

The reported co-expression analysis identified modules associated with these lncRNAs and highlighted biological processes relevant to metastatic CRPC, including:

- epithelial–mesenchymal transition;
- mesenchymal–epithelial transition;
- purine metabolism;
- transcriptional regulation;
- immune-system processes.

The study also reported the following hub-gene candidates:

- `SOD2`
- `PRKCA`
- `IL6`
- `ITGAM`

## Aim of This Repository

This repository provides a structured computational workflow to:

- download the `GSE74685` dataset from the Gene Expression Omnibus;
- extract the expression matrix and sample metadata;
- classify samples according to metastatic site;
- filter low-expression probes;
- compare Bone and Visceral metastasis groups using `limma`;
- construct a signed WGCNA co-expression network;
- identify gene co-expression modules;
- calculate module–trait correlations;
- calculate module-membership values;
- combine differential-expression results with WGCNA results;
- identify candidate lncRNAs using compatible platform annotation.

## Important Scope Note

The current scripts perform a Bone-versus-Visceral comparison based on metastatic-site information inferred from GEO sample titles.

This workflow is an adaptation inspired by the reference article. It should not be described as an exact reproduction of every analysis reported in the paper.

The automatic sample-classification procedure searches the sample-title field for terms such as:

- `bone`
- `visceral`
- `liver`
- `lymph`

Metastatic-site classifications must be manually verified against the original GEO metadata before interpreting the results.

## Repository Structure

```text
CRPC-lncRNA-WGCNA-Analysis-of-GSE74685/
├── README.md
├── LICENSE
├── scripts/
│   ├── Script 01_load_packages.R
│   ├── Script 02_download_preprocess_GSE74685.R
│   ├── Script 03_DE_Bone_vs_Visceral.R
│   ├── Script 04_WGCNA_modules_traits.R
│   └── Script 05_extract_lncRNA_candidates.R
├── data/
├── results/
├── figures/
└── docs/
```

## Requirements

The following resources are recommended:

- R version 4.2 or newer;
- internet access for the initial package installation and GEO download;
- sufficient memory for WGCNA network construction;
- a compatible probe-annotation file for identifying lncRNAs.

## Required R Packages

### Bioconductor Packages

- `GEOquery`
- `Biobase`
- `limma`

### CRAN Packages

- `WGCNA`
- `dynamicTreeCut`
- `tidyverse`

Missing dependencies are installed automatically by `Script 01_load_packages.R`.

## Running the Workflow

Start R or RStudio from the root directory of the repository.

Do not set the `scripts` folder as the working directory.

Run the scripts in numerical order:

```r
source("scripts/Script 01_load_packages.R")
source("scripts/Script 02_download_preprocess_GSE74685.R")
source("scripts/Script 03_DE_Bone_vs_Visceral.R")
source("scripts/Script 04_WGCNA_modules_traits.R")
source("scripts/Script 05_extract_lncRNA_candidates.R")
```

## Pipeline Stages

### Stage 01: Package Setup

The `Script 01_load_packages.R` script:

- installs missing CRAN packages;
- installs missing Bioconductor packages;
- loads the required packages;
- disables automatic string-to-factor conversion;
- enables WGCNA multithreading when supported.

### Stage 02: GEO Download and Preprocessing

The `Script 02_download_preprocess_GSE74685.R` script:

- downloads `GSE74685` using `GEOquery`;
- extracts the expression matrix;
- extracts the phenotype and sample-metadata table;
- aligns phenotype rows with expression-matrix columns;
- derives the `site_group` variable;
- classifies samples as Bone, Visceral, LymphNode, or Other;
- filters low-expression probes;
- saves the processed data in the `data` directory.

The default filtering rule retains probes with expression greater than `5` in at least 20% of samples.

### Stage 03: Differential-Expression Analysis

The `Script 03_DE_Bone_vs_Visceral.R` script:

- retains Bone and Visceral samples;
- aligns expression data with phenotype information;
- constructs the design matrix;
- defines the `Bone - Visceral` contrast;
- fits a linear model using `limma`;
- applies empirical Bayes moderation;
- generates a ranked differential-expression table;
- saves the results in the `results` directory.

A positive `logFC` indicates higher expression in Bone samples relative to Visceral samples.

A negative `logFC` indicates higher expression in Visceral samples relative to Bone samples.

### Stage 04: WGCNA Module and Trait Analysis

The `Script 04_WGCNA_modules_traits.R` script:

- transposes the expression matrix to samples × probes;
- identifies unsuitable samples and probes using `goodSamplesGenes`;
- evaluates candidate soft-thresholding powers;
- selects an appropriate soft-thresholding power;
- constructs a signed co-expression network;
- detects co-expression modules;
- merges similar modules;
- calculates module eigengenes;
- creates Bone and Visceral trait variables;
- calculates module–trait correlations;
- saves WGCNA objects and correlation results.

### Stage 05: Candidate Gene and lncRNA Extraction

The `Script 05_extract_lncRNA_candidates.R` script:

- loads the WGCNA objects;
- loads the differential-expression results;
- calculates module-membership values;
- combines WGCNA and differential-expression statistics;
- optionally adds gene symbols and transcript biotypes;
- identifies Bone-associated modules;
- applies module-membership and differential-expression filters;
- extracts candidate probes;
- extracts candidate lncRNAs;
- checks the key lncRNAs reported in the reference article.
