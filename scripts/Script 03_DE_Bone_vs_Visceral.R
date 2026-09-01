# Project: CRPC lncRNA WGCNA (GSE74685)
# Goal: Differential-expression analysis of Bone versus Visceral samples.

source(file.path("scripts", "Script 01_load_packages.R"))

gse_id <- "GSE74685"
expr_file <- file.path("data", paste0(gse_id, "_expr_filtered.rds"))
pheno_file <- file.path("data", paste0(gse_id, "_pheno.csv"))
deg_outfile <- file.path("results", paste0(
  gse_id, "_DEG_Bone_vs_Visceral_limma.csv"))

if (!file.exists(expr_file)) stop("Run Script 02 first: missing ", expr_file)
if (!file.exists(pheno_file)) stop("Run Script 02 first: missing ", pheno_file)

expr_filt <- readRDS(expr_file)
pheno_df <- read.csv(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample_id", "site_group") %in% names(pheno_df))) {
  stop("Phenotype table must contain sample_id and site_group.")
}
rownames(pheno_df) <- pheno_df$sample_id

keep_samples <- tolower(pheno_df$site_group) %in% c("bone", "visceral")
pheno_sub <- pheno_df[keep_samples, , drop = FALSE]
common_samples <- colnames(expr_filt)[colnames(expr_filt) %in%
                                      rownames(pheno_sub)]
if (!length(common_samples)) stop("Expression and phenotype IDs do not overlap.")

pheno_sub <- pheno_sub[common_samples, , drop = FALSE]
expr_sub <- expr_filt[, common_samples, drop = FALSE]
group <- factor(pheno_sub$site_group, levels = c("Visceral", "Bone"))
if (anyNA(group) || any(table(group) < 2)) {
  stop("Bone and Visceral must each contain at least two samples.")
}
print(table(group))

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contrast_matrix <- limma::makeContrasts(
  Bone_vs_Visceral = Bone - Visceral, levels = design)
fit <- limma::lmFit(expr_sub, design)
fit <- limma::contrasts.fit(fit, contrast_matrix)
fit <- limma::eBayes(fit)
deg_table <- limma::topTable(fit, coef = "Bone_vs_Visceral",
                             number = Inf, sort.by = "P")

dir.create("results", recursive = TRUE, showWarnings = FALSE)
write.csv(deg_table, deg_outfile, row.names = TRUE)
message("Saved: ", deg_outfile)

