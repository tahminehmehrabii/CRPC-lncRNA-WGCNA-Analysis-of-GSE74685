# Project: CRPC lncRNA WGCNA (GSE74685)
# Goal: Build WGCNA modules and calculate module-trait correlations.

source(file.path("scripts", "Script 01_load_packages.R"))

gse_id <- "GSE74685"
expr_file <- file.path("data", paste0(gse_id, "_expr_filtered.rds"))
pheno_file <- file.path("data", paste0(gse_id, "_pheno.csv"))
if (!file.exists(expr_file) || !file.exists(pheno_file)) {
  stop("Run Script 02 before Script 04.")
}

expr_filt <- readRDS(expr_file)
pheno_df <- read.csv(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)
rownames(pheno_df) <- pheno_df$sample_id
common_samples <- colnames(expr_filt)[colnames(expr_filt) %in%
                                      rownames(pheno_df)]
datExpr <- as.data.frame(t(expr_filt[, common_samples, drop = FALSE]))
pheno_df <- pheno_df[common_samples, , drop = FALSE]

gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  pheno_df <- pheno_df[rownames(datExpr), , drop = FALSE]
}

powers <- c(1:10, seq(12, 30, 2))
sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = powers,
                                networkType = "signed", verbose = 5)
fit <- sft$fitIndices
eligible <- fit$Power[is.finite(fit$SFT.R.sq) & fit$SFT.R.sq >= 0.80]
soft_power <- if (length(eligible)) eligible[1] else
  fit$Power[which.max(replace(fit$SFT.R.sq, !is.finite(fit$SFT.R.sq), -Inf))]
if (!length(soft_power) || !is.finite(soft_power)) soft_power <- 6
message("Selected soft-threshold power: ", soft_power)

net <- WGCNA::blockwiseModules(
  datExpr, power = soft_power, TOMType = "signed",
  minModuleSize = 30, mergeCutHeight = 0.25,
  numericLabels = FALSE, pamRespectsDendro = FALSE,
  saveTOMs = FALSE, verbose = 3
)
moduleColors <- net$colors
names(moduleColors) <- colnames(datExpr)
MEs <- WGCNA::orderMEs(net$MEs)

traitData <- data.frame(
  Bone = as.numeric(pheno_df$site_group == "Bone"),
  Visceral = as.numeric(pheno_df$site_group == "Visceral"),
  row.names = rownames(pheno_df)
)
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitP <- WGCNA::corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))
traitCor_df <- data.frame(
  Module = sub("^ME", "", rownames(moduleTraitCor)),
  Cor_Bone = moduleTraitCor[, "Bone"],
  P_Bone = moduleTraitP[, "Bone"],
  Cor_Visceral = moduleTraitCor[, "Visceral"],
  P_Visceral = moduleTraitP[, "Visceral"],
  row.names = NULL
)

dir.create("data", recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)
saveRDS(datExpr, file.path("data", paste0(gse_id, "_WGCNA_datExpr.rds")))
saveRDS(moduleColors, file.path("data", paste0(gse_id,
                                                "_WGCNA_moduleColors.rds")))
saveRDS(MEs, file.path("data", paste0(gse_id, "_WGCNA_MEs.rds")))
write.csv(traitCor_df, file.path("results", paste0(
  gse_id, "_WGCNA_module_trait_correlations.csv")), row.names = FALSE)
message("Finished Script 04.")

