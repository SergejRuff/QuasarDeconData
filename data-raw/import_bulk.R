rm(list=ls())


bulk_path <- "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/real_bulk_pbmc_8k/bulk/rds_files"


GSE107011 <- readRDS(list.files(bulk_path,full.names = TRUE)[[1]])

GSE107572 <- readRDS(list.files(bulk_path,full.names = TRUE)[[2]])


GSE120502 <- readRDS(list.files(bulk_path,full.names = TRUE)[[3]])


GSE65133 <- readRDS(list.files(bulk_path,full.names = TRUE)[[4]])

library(illuminaHumanv4.db)
symbols <- mapIds(illuminaHumanv4.db,
                  keys = rownames(GSE65133$bulk_expression_profiles),
                  column = "SYMBOL", keytype = "PROBEID",
                  multiVals = "first")

X <- GSE65133$bulk_expression_profiles

keep <- !is.na(symbols)
X <- X[keep, , drop = FALSE]
sym <- symbols[keep]

v <- apply(X, 1, var)
best <- tapply(seq_along(sym), sym, function(i) i[which.max(v[i])])

X <- X[unlist(best), , drop = FALSE]
rownames(X) <- names(best)


GSE65133$bulk_expression_profiles <- X

usethis::use_data(
  GSE107011,
  overwrite = TRUE,
  compress = "xz"
)

usethis::use_data(
  GSE107572,
  overwrite = TRUE,
  compress = "xz"
)

usethis::use_data(
  GSE120502,
  overwrite = TRUE,
  compress = "xz"
)

usethis::use_data(
  GSE65133,
  overwrite = TRUE,
  compress = "xz"
)

tools::resaveRdaFiles("data/", compress = "auto")
tools::checkRdaFiles("data/")  # shows which compression won per file


names(GSE65133)

dim(GSE65133$bulk_expression_profiles)

colnames(GSE65133$ground_truth_proportions)
dim(GSE65133$ground_truth_proportions)



####################################################

rm(list = ls())

train_pb <- "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/real_bulk_pbmc_8k/pseudo_bulk_train/pbmc_healthy_pb.RDS"

keep_names <- c("bulk_expression_profiles", "ground_truth_proportions",
                "cells_per_sample")

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

train_pbmc <- readRDS(train_pb)[keep_names]

expr <- train_pbmc$bulk_expression_profiles
prop <- train_pbmc$ground_truth_proportions
cps  <- train_pbmc$cells_per_sample

n_samples <- nrow(prop)
stopifnot(ncol(expr) == n_samples || nrow(expr) == n_samples)
samples_in_cols <- ncol(expr) == n_samples

cat("expr:", nrow(expr), "x", ncol(expr), "| samples in",
    if (samples_in_cols) "cols" else "rows", "\n")
cat("prop:", nrow(prop), "x", ncol(prop), "\n")
cat("cell types:", paste(colnames(prop), collapse = ", "), "\n")

n_chunks <- 10L
idx_list <- split(seq_len(n_samples), cut(seq_len(n_samples), n_chunks, labels = FALSE))

for (i in seq_along(idx_list)) {
  idx <- idx_list[[i]]

  chunk <- list(
    bulk_expression_profiles = if (samples_in_cols) {
      expr[, idx, drop = FALSE]
    } else {
      expr[idx, , drop = FALSE]
    },
    ground_truth_proportions = prop[idx, , drop = FALSE],
    cells_per_sample = if (length(cps) == n_samples) cps[idx] else cps
  )
  attr(chunk, "sample_axis") <- if (samples_in_cols) "cols" else "rows"

  f <- file.path("inst/extdata", sprintf("pbmc_healthy_pb_%d.rds", i))
  con <- xzfile(f, compression = 9)
  saveRDS(chunk, con)
  close(con)

  cat(sprintf("%s: %.1f MB\n", basename(f), file.info(f)$size / 1024^2))
}