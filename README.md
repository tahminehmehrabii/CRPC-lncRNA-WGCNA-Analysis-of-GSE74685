# CRPC-lncRNA-WGCNA-Analysis-of-GSE74685

CRPC lncRNA WGCNA Analysis of GSE74685

This repository contains an R-based workflow for differential-expression and weighted gene co-expression network analysis of the GSE74685 dataset. The workflow was developed in connection with the study of dysregulated long non-coding RNAs and epithelial–mesenchymal transition in castration-resistant prostate cancer.

Reference article

Mehrabi T, Heidarzadehpilehrood R, Mobasheri M, Sobati T, Heshmati M, Pirhoushiaran M. Dysregulated key long non-coding RNAs TP53TG1, RFPL1S, DLEU1, and HCG4 associated with epithelial-mesenchymal transition (EMT) in castration-resistant prostate cancer. Advances in Cancer Biology - Metastasis. 2025;13:100132.

DOI: 10.1016/j.adcanc.2025.100132

Scientific background

Castration-resistant prostate cancer is an advanced form of prostate cancer with a high risk of progression and metastasis. Long non-coding RNAs can contribute to cancer progression by affecting gene regulation, cellular plasticity, invasion, and epithelial–mesenchymal transition.

The reference article analyzed GSE74685 and reported four dysregulated lncRNAs:

TP53TG1

RFPL1S

DLEU1

HCG4

The reported co-expression analysis identified modules associated with these lncRNAs and highlighted pathways relevant to metastasis, including epithelial/mesenchymal transition, purine metabolism, transcriptional regulation, and immune-system processes. The article also reported SOD2, PRKCA, IL6, and ITGAM as hub-gene candidates.

Aim of this repository

The pipeline provides a structured starting point to:

download and preprocess the GSE74685 expression dataset;

classify samples by metastatic site using GEO sample titles;

compare Bone and Visceral metastasis groups with limma;

construct signed WGCNA co-expression modules;

calculate module–trait correlations;

combine differential-expression statistics with module membership;

extract candidate lncRNAs when compatible probe annotation is available.

Important scope note

The present scripts perform a Bone-versus-Visceral comparison inferred from the GEO sample titles. This is an adaptation inspired by the reference article and should not be described as an exact reproduction of every analysis in the paper.

Metastatic-site labels should be manually checked before interpreting the results. The automatic grouping rules search the sample-title field for terms such as bone, visceral, liver, and lymph.

Repository structure

CRPC_lncRNA_WGCNA/
├── README.md
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

Requirements

R version 4.2 or newer is recommended.

Internet access is required for the first GEO download and package installation.

Sufficient memory is required for WGCNA network construction.

The workflow uses the following packages:

Bioconductor: GEOquery, Biobase, and limma

CRAN: WGCNA, dynamicTreeCut, and tidyverse

Missing dependencies are installed automatically by Script 01.

Running the workflow

Start R from the repository root, not from the scripts directory. Run the scripts in numerical order:

source("scripts/Script 01_load_packages.R")
source("scripts/Script 02_download_preprocess_GSE74685.R")
source("scripts/Script 03_DE_Bone_vs_Visceral.R")
source("scripts/Script 04_WGCNA_modules_traits.R")
source("scripts/Script 05_extract_lncRNA_candidates.R")

Pipeline stages

Stage 01: package setup

Script 01_load_packages.R installs missing CRAN and Bioconductor packages, loads them, disables automatic string-to-factor conversion, and enables WGCNA multithreading where supported.

Stage 02: download and preprocessing

Script 02_download_preprocess_GSE74685.R:

downloads GSE74685 with GEOquery;

extracts the expression matrix and phenotype table;

aligns phenotype rows with expression columns;

derives the site_group phenotype;

retains probes with expression greater than 5 in at least 20% of samples;

saves cached inputs in data/.

Stage 03: differential expression

Script 03_DE_Bone_vs_Visceral.R keeps Bone and Visceral samples, fits a limma model, evaluates the Bone - Visceral contrast, and saves a ranked table containing log-fold changes, P values, and multiple-testing-adjusted P values.

A positive logFC indicates higher expression in Bone samples relative to Visceral samples.

Stage 04: WGCNA

Script 04_WGCNA_modules_traits.R:

transposes the expression matrix to samples × probes;

removes unsuitable samples or probes using goodSamplesGenes;

evaluates candidate soft-thresholding powers;

constructs a signed co-expression network;

detects and merges co-expression modules;

calculates module eigengenes;

correlates module eigengenes with Bone and Visceral traits.

Stage 05: candidate extraction

Script 05_extract_lncRNA_candidates.R combines module assignments, module membership values, differential-expression results, and optional platform annotation. Candidate probes are selected using module–trait association, module membership, log-fold-change, and adjusted-P-value criteria.
