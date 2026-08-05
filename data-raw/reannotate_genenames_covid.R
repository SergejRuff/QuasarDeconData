rm(list = ls())

library(Seurat)
library(Matrix)
library(org.Hs.eg.db)
library(AnnotationDbi)

split_dir  <- "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/covid_immune_atlas/split/"
pb_all_dir <- "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/covid_immune_atlas/pseudo_bulk_all"
pb_train_f <- "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/covid_immune_atlas/pseudo_bulk_train/train_bulk.rds"

train_f <- file.path(split_dir, "train_cov_imm.RDS")
test_f  <- file.path(split_dir, "test_cov_imm.RDS")
pb_files <- list.files(pb_all_dir, full.names = TRUE, pattern = "\\.rds$")

sym <- unname(mapIds(org.Hs.eg.db,
                     keys = sub("\\..*$", "", rownames(readRDS(train_f))),
                     keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
keep <- !is.na(sym) & sym != ""
grp  <- fac2sparse(factor(sym[keep]))

collapse <- function(m) grp %*% m[keep, , drop = FALSE]

convert_seurat <- function(f, proj) {
  o <- readRDS(f)
  meta <- o[[]]
  n <- CreateSeuratObject(counts = collapse(LayerData(o, assay = "RNA", layer = "counts")),
                          meta.data = meta, assay = "RNA", project = proj,
                          min.cells = 0, min.features = 0)
  n <- NormalizeData(n, verbose = FALSE)
  for (r in Reductions(o)) n[[r]] <- o[[r]]
  if ("cell_type" %in% colnames(meta)) Idents(n) <- n$cell_type
  saveRDS(n, f)
  rm(o, n); gc()
}

convert_bulk <- function(f) {
  x <- readRDS(f)
  x$bulk_expression_profiles <- as.matrix(collapse(as.matrix(x$bulk_expression_profiles)))
  saveRDS(x, f)
  rm(x); gc()
}

convert_seurat(train_f, "train_cov")
convert_seurat(test_f, "test_cov")
for (f in pb_files) convert_bulk(f)
convert_bulk(pb_train_f)