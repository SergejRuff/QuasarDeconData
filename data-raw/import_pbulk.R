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

# keep the same two slots as the test data (drop if not wanted)
train_pb <- train_pb[keep_names]

format(object.size(train_pb), units = "auto")

con <- xzfile(
  file.path("inst/extdata", "training_pb_cov.rds"),
  compression = 9
)
saveRDS(train_pb, con)
close(con)