rm(list = ls())

library(Seurat)

base_folder <- "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/covid_immune_atlas/split/"

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

prep <- function(file, project) {
  obj <- readRDS(file)

  cnt  <- LayerData(obj, assay = "RNA", layer = "counts")
  meta <- obj[[]]

  out <- CreateSeuratObject(
    counts    = cnt,
    meta.data = meta,
    assay     = "RNA",
    project   = project
  )

  out@meta.data <- meta[, c("cell_type", "donor_id"), drop = FALSE]
  out$cell_type <- droplevels(factor(as.character(out$cell_type)))
  Idents(out) <- out$cell_type

  rm(obj, cnt); gc()
  out
}

pack <- function(obj, file) {
  m <- LayerData(obj, assay = "RNA", layer = "counts")
  m <- as(m, "dgCMatrix")
  payload <- list(
    i = m@i, p = m@p,
    x = writeBin(m@x, raw(), size = 4L),
    n = length(m@x),
    Dim = m@Dim, dn = m@Dimnames,
    meta = obj@meta.data
  )
  con <- xzfile(file, compression = 9); saveRDS(payload, con); close(con)
  file.info(file)$size / 1024^2
}

train_cov <- prep(paste0(base_folder, "train_cov_imm.RDS"), "train_cov")
pack(train_cov, "inst/extdata/train_cov.rds")
rm(train_cov); gc()

test_cov <- prep(paste0(base_folder, "test_cov_imm.RDS"), "test_cov")
pack(test_cov, "inst/extdata/test_cov.rds")
rm(test_cov); gc()

healthy_cov <- prep(paste0(base_folder, "train_healthy_ref.RDS"), "healthy_cov")
print(table(healthy_cov$cell_type))
pack(healthy_cov, "inst/extdata/healthy_cov.rds")
rm(healthy_cov); gc()

for (f in c("train_cov.rds", "test_cov.rds", "healthy_cov.rds")) {
  cat(f, ":", round(file.info(file.path("inst/extdata", f))$size / 1024^2, 1), "MB\n")
}