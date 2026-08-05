rm(list=ls())

library(Seurat)

base_folder = "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/covid_immune_atlas/split/"

train_cov <- paste0(base_folder,"train_cov_imm.RDS")
test_cov <- paste0(base_folder,"test_cov_imm.RDS")

train_cov <- readRDS(train_cov)
test_cov <- readRDS(test_cov)



train_cov_counts <- LayerData(
  train_cov,
  assay = "RNA",
  layer = "counts"
)

test_cov_counts <- LayerData(
  test_cov,
  assay = "RNA",
  layer = "counts"
)

# Extract metadata
train_meta_data <- train_cov[[]]
test_meta_data <- test_cov[[]]

# Rebuild object with only counts
train_cov <- CreateSeuratObject(
  counts = train_cov_counts,
  meta.data = train_meta_data,
  assay = "RNA",
  project = "train_cov"
)

test_cov <- CreateSeuratObject(
  counts = test_cov_counts,
  meta.data = test_meta_data,
  assay = "RNA",
  project = "test_cov"
)

train_cov@meta.data <- train_meta_data[, c("cell_type", "donor_id"), drop = FALSE]
test_cov@meta.data <- test_meta_data[, c("cell_type", "donor_id"), drop = FALSE]

# Optional: preserve cell identities if you used celltype as identities
Idents(train_cov) <- train_cov$cell_type
Idents(test_cov) <- test_cov$cell_type




dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

# con <- xzfile("inst/extdata/train_cov.rds", compression = 9)
# saveRDS(train_cov, con); close(con)

# con <- xzfile("inst/extdata/test_cov.rds", compression = 9)
# saveRDS(test_cov, con); close(con)


pack <- function(obj, file) {
  m <- LayerData(obj, assay = "RNA", layer = "counts")
  payload <- list(
    i = m@i, p = m@p,
    x = writeBin(m@x, raw(), size = 4L),   # float64 -> float32 bytes
    n = length(m@x),
    Dim = m@Dim, dn = m@Dimnames,
    meta = obj@meta.data
  )
  con <- xzfile(file, compression = 9); saveRDS(payload, con); close(con)
  file.info(file)$size / 1024^2
}
pack(train_cov, "inst/extdata/train_cov.rds")
pack(test_cov,  "inst/extdata/test_cov.rds")
# saveRDS(
#   train_cov,
#   file = "inst/extdata/train_cov.rds",
#   compress = "xz"
# )

# Check compressed size
print(file.info("inst/extdata/train_cov.rds")$size / 1024^2)


# saveRDS(
#   test_cov,
#   file = "inst/extdata/test_cov.rds",
#   compress = "xz"
# )

# Check compressed size
print(file.info("inst/extdata/test_cov.rds")$size / 1024^2)
