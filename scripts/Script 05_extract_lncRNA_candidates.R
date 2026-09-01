# Project: CRPC lncRNA WGCNA (GSE74685)
# Goal: Combine WGCNA membership, DE results, and optional annotation.

source(file.path("scripts", "Script 01_load_packages.R"))
gse_id <- "GSE74685"

datExpr_file <- file.path("data", paste0(gse_id, "_WGCNA_datExpr.rds"))
moduleColors_file <- file.path("data", paste0(gse_id,
                                               "_WGCNA_moduleColors.rds"))
MEs_file <- file.path("data", paste0(gse_id, "_WGCNA_MEs.rds"))
traitCor_file <- file.path("results", paste0(
  gse_id, "_WGCNA_module_trait_correlations.csv"))
deg_file <- file.path("results", paste0(
  gse_id, "_DEG_Bone_vs_Visceral_limma.csv"))
required <- c(datExpr_file, moduleColors_file, MEs_file, traitCor_file, deg_file)
if (any(!file.exists(required))) stop("Missing inputs; run Scripts 02–04 first.")

datExpr <- readRDS(datExpr_file)
moduleColors <- readRDS(moduleColors_file)
MEs <- readRDS(MEs_file)
traitCor_df <- read.csv(traitCor_file, stringsAsFactors = FALSE)
deg_tab <- read.csv(deg_file, row.names = 1, check.names = FALSE)
if (!identical(colnames(datExpr), names(moduleColors))) {
  moduleColors <- moduleColors[colnames(datExpr)]
}

kME_mat <- as.data.frame(cor(datExpr, MEs, use = "p"))
colnames(kME_mat) <- paste0("kME_", colnames(kME_mat))
geneInfo <- data.frame(Probe = colnames(datExpr), Module = unname(moduleColors),
                       kME_mat, stringsAsFactors = FALSE, check.names = FALSE)
rownames(geneInfo) <- geneInfo$Probe
common_ids <- intersect(rownames(deg_tab), rownames(geneInfo))
if (!length(common_ids)) stop("DEG and WGCNA probe IDs do not overlap.")
geneInfo_full <- cbind(geneInfo[common_ids, , drop = FALSE],
                       deg_tab[common_ids, , drop = FALSE])

annot_file <- file.path("docs", "GPL_annotation.csv")
if (file.exists(annot_file)) {
  annot <- read.csv(annot_file, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("ProbeID", "GeneSymbol", "Biotype")
  if (!all(needed %in% names(annot))) stop("Annotation needs: ",
                                           paste(needed, collapse = ", "))
  annot <- annot[, needed]
  names(annot)[1] <- "Probe"
  geneInfo_full <- merge(geneInfo_full, annot, by = "Probe", all.x = TRUE,
                         sort = FALSE)
} else {
  message("No docs/GPL_annotation.csv; lncRNA-specific output will be empty.")
}

cor_cut <- 0.30
p_cut <- 0.05
kME_cut <- 0.5
logFC_cut <- 1
fdr_cut <- 0.05
sig_modules <- traitCor_df$Module[
  abs(traitCor_df$Cor_Bone) >= cor_cut & traitCor_df$P_Bone <= p_cut]
kme_cols <- grep("^kME_", names(geneInfo_full), value = TRUE)
p_col <- if ("adj.P.Val" %in% names(geneInfo_full)) "adj.P.Val" else "P.Value"
own_kme <- vapply(seq_len(nrow(geneInfo_full)), function(i) {
  col <- paste0("kME_ME", geneInfo_full$Module[i])
  if (!col %in% kme_cols) return(FALSE)
  is.finite(geneInfo_full[[col]][i]) && abs(geneInfo_full[[col]][i]) >= kME_cut
}, logical(1))
keep <- geneInfo_full$Module %in% sig_modules & own_kme &
  abs(geneInfo_full$logFC) >= logFC_cut &
  geneInfo_full[[p_col]] <= fdr_cut
gene_candidates <- geneInfo_full[keep, , drop = FALSE]

if ("Biotype" %in% names(gene_candidates)) {
  lnc_candidates <- gene_candidates[
    grepl("lnc|long non.?coding", gene_candidates$Biotype,
          ignore.case = TRUE), , drop = FALSE]
} else lnc_candidates <- gene_candidates[0, , drop = FALSE]

key_symbols <- c("TP53TG1", "RFPL1S", "DLEU1")
key_hits <- if ("GeneSymbol" %in% names(geneInfo_full)) {
  geneInfo_full[geneInfo_full$GeneSymbol %in% key_symbols, , drop = FALSE]
} else geneInfo_full[0, , drop = FALSE]

dir.create("results", recursive = TRUE, showWarnings = FALSE)
write.csv(geneInfo_full, file.path("results", paste0(
  gse_id, "_WGCNA_geneInfo_allProbes.csv")), row.names = FALSE)
write.csv(gene_candidates, file.path("results", paste0(
  gse_id, "_WGCNA_moduleCandidates.csv")), row.names = FALSE)
write.csv(lnc_candidates, file.path("results", paste0(
  gse_id, "_WGCNA_lncRNA_candidates.csv")), row.names = FALSE)
write.csv(key_hits, file.path("results", paste0(
  gse_id, "_WGCNA_key_lncRNAs_from_paper.csv")), row.names = FALSE)
message("Finished Script 05; candidates: ", nrow(gene_candidates),
        ", lncRNAs: ", nrow(lnc_candidates), ".")

