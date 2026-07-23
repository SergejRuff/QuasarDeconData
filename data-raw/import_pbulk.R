rm(list=ls())

library(Seurat)

base_folder = "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/covid_immune_atlas/pseudo_bulk_all"


rds_files <- list.files(
  base_folder,
  full.names = TRUE,
  pattern = "\\.rds$"
)

pbulk_list <- lapply(rds_files, readRDS)


keep_names <- c("bulk_expression_profiles", "ground_truth_proportions")

pbulk_list <- lapply(pbulk_list, function(x) {
  x[keep_names]
})

format(object.size(pbulk_list), units = "auto")



dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(pbulk_list)) {
  con <- xzfile(
    file.path("inst/extdata", paste0("cov_pbulk_", i, ".rds")),
    compression = 9
  )
  
  saveRDS(pbulk_list[[i]], con)
  close(con)
}


train_file <- "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/covid_immune_atlas/pseudo_bulk_train/train_bulk.rds"

train_pb <- readRDS(train_file)
keep_names <- c("bulk_expression_profiles", "ground_truth_proportions")
train_pb <- train_pb[keep_names]

expr <- train_pb$bulk_expression_profiles
prop <- train_pb$ground_truth_proportions

n_samples <- nrow(prop)
stopifnot(ncol(expr) == n_samples || nrow(expr) == n_samples)
samples_in_cols <- ncol(expr) == n_samples

n_chunks <- 10L
idx_list <- split(seq_len(n_samples), cut(seq_len(n_samples), n_chunks, labels = FALSE))

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

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

  con <- xzfile(
    file.path("inst/extdata", sprintf("training_pb_cov_%d.rds", i)),
    compression = 9
  )
  saveRDS(chunk, con)
  close(con)
}

# sanity check — every file must be well under 100 MB
sizes <- file.size(list.files("inst/extdata", pattern = "^training_pb_cov_", full.names = TRUE))
round(sizes / 1024^2, 1)