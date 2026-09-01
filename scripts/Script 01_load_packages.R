# Project: CRPC lncRNA WGCNA (GSE74685)
# Goal: Install and load required packages.

bioc_packages <- c("GEOquery", "Biobase", "limma")
cran_packages <- c("WGCNA", "dynamicTreeCut", "tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_bioc <- bioc_packages[!vapply(
  bioc_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

missing_cran <- cran_packages[!vapply(
  cran_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

invisible(lapply(c(bioc_packages, cran_packages), library,
                 character.only = TRUE))

options(stringsAsFactors = FALSE)
WGCNA::allowWGCNAThreads()
message("All required packages are loaded.")

