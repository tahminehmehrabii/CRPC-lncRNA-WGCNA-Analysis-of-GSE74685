# Project: CRPC lncRNA WGCNA (GSE74685)
# Goal: Download, harmonise, filter, and save expression/phenotype data.

source(file.path("scripts", "Script 01_load_packages.R"))

dirs <- c("data", "results", "figures", "docs")
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

gse_id <- "GSE74685"
expr_outfile_raw  <- file.path("data", paste0(gse_id, "_expr_raw.rds"))
expr_outfile_filt <- file.path("data", paste0(gse_id, "_expr_filtered.rds"))
pheno_outfile     <- file.path("data", paste0(gse_id, "_pheno.csv"))
expr_threshold <- 5

needed <- c(expr_outfile_raw, expr_outfile_filt, pheno_outfile)
if (!all(file.exists(needed))) {
  message("Downloading ", gse_id, " from GEO...")
  gse_list <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE)
  if (!length(gse_list)) stop("GEO returned no ExpressionSet for ", gse_id)
  gset <- gse_list[[1]]

  expr_mat <- Biobase::exprs(gset)
  pheno_df <- Biobase::pData(gset)
  sample_ids <- colnames(expr_mat)

  pheno_df$sample_id <- rownames(pheno_df)
  match_index <- match(sample_ids, pheno_df$sample_id)
  if (anyNA(match_index)) stop("Some expression samples lack phenotype rows.")
  pheno_df <- pheno_df[match_index, , drop = FALSE]
  rownames(pheno_df) <- pheno_df$sample_id

  title_vec <- as.character(pheno_df$title)
  site_group <- dplyr::case_when(
    grepl("bone", title_vec, ignore.case = TRUE) ~ "Bone",
    grepl("viscer|liver", title_vec, ignore.case = TRUE) ~ "Visceral",
    grepl("lymph", title_vec, ignore.case = TRUE) ~ "LymphNode",
    TRUE ~ "Other"
  )
  pheno_df$site_group <- factor(site_group)

  keep_probes <- rowMeans(expr_mat > expr_threshold, na.rm = TRUE) >= 0.20
  expr_filt <- expr_mat[keep_probes, , drop = FALSE]

  saveRDS(expr_mat, expr_outfile_raw)
  saveRDS(expr_filt, expr_outfile_filt)
  write.csv(pheno_df, pheno_outfile, row.names = FALSE)
} else {
  message("Cached files found; skipping download.")
  expr_filt <- readRDS(expr_outfile_filt)
  pheno_df <- read.csv(pheno_outfile, stringsAsFactors = FALSE,
                       check.names = FALSE)
  rownames(pheno_df) <- pheno_df$sample_id
}

message("Filtered matrix: ", nrow(expr_filt), " probes x ",
        ncol(expr_filt), " samples.")
message("Finished Script 02.")

