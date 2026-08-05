rm(list = ls())

library(Seurat)

base_folder <- "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/covid_immune_atlas/pseudo_bulk_all"
train_file   <- "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/covid_immune_atlas/pseudo_bulk_train/train_bulk.rds"
healthy_file <- "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/covid_immune_atlas/split/train_healthy_pb.RDS"

keep_names <- c("bulk_expression_profiles", "ground_truth_proportions")

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(base_folder, full.names = TRUE, pattern = "\\.rds$")
pbulk_list <- lapply(rds_files, function(f) readRDS(f)[keep_names])

format(object.size(pbulk_list), units = "auto")

for (i in seq_along(pbulk_list)) {
  con <- xzfile(file.path("inst/extdata", paste0("cov_pbulk_", i, ".rds")),
                compression = 9)
  saveRDS(pbulk_list[[i]], con)
  close(con)
}

rm(pbulk_list); gc()

chunk_and_save <- function(file, prefix, n_chunks = 10L) {
  pb <- readRDS(file)[keep_names]

  expr <- pb$bulk_expression_profiles
  prop <- pb$ground_truth_proportions

  n_samples <- nrow(prop)
  stopifnot(ncol(expr) == n_samples || nrow(expr) == n_samples)
  samples_in_cols <- ncol(expr) == n_samples

  idx_list <- split(seq_len(n_samples), cut(seq_len(n_samples), n_chunks, labels = FALSE))

  for (i in seq_along(idx_list)) {
    idx <- idx_list[[i]]

    chunk <- list(
      bulk_expression_profiles = if (samples_in_cols) {
        expr[, idx, drop = FALSE]
      } else {
        expr[idx, , drop = FALSE]
      },
      ground_truth_proportions = prop[idx, , drop = FALSE]
    )
    attr(chunk, "sample_axis") <- if (samples_in_cols) "cols" else "rows"

    con <- xzfile(file.path("inst/extdata", sprintf("%s_%d.rds", prefix, i)),
                  compression = 9)
    saveRDS(chunk, con)
    close(con)
  }

  rm(pb, expr, prop); gc()
  invisible(n_samples)
}

chunk_and_save(train_file,   "training_pb_cov")
chunk_and_save(healthy_file, "healthy_pb_cov")

sizes <- file.size(list.files("inst/extdata", pattern = "\\.rds$", full.names = TRUE))
data.frame(file = basename(list.files("inst/extdata", pattern = "\\.rds$")),
           mb = round(sizes / 1024^2, 1))