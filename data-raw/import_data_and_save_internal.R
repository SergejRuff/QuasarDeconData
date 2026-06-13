rm(list=ls())

rds_file_ref <- "/home/sergej/Schreibtisch/program/10_compare_deconvolution_methods/analysis/data_derived/real_bulk_pbmc_8k/ref/processed_seurat/pbmc6k_8k_10k_combined_processed.rds"

pbmc_sc_ref <- readRDS(rds_file_ref)

pbmc_sc_ref <- subset(pbmc_sc_ref, subset = celltype != "Dendritics")
pbmc_sc_ref$celltype <- droplevels(pbmc_sc_ref$celltype)


counts_mat <- pbmc_sc_ref[["RNA"]]$counts

# Keep genes expressed in at least one cell
# genes_keep <- rownames(counts_mat)[Matrix::rowSums(counts_mat > 0) > 0]

# length(genes_keep)

# pbmc_sc_ref <- subset(
#   pbmc_sc_ref,
#   features = genes_keep
# )

table(pbmc_sc_ref$donor_id)
table(pbmc_sc_ref$celltype)
pbmc_sc_ref


library(Seurat)


counts_mat <- LayerData(
  pbmc_sc_ref,
  assay = "RNA",
  layer = "counts"
)

# Extract metadata
meta_data <- pbmc_sc_ref[[]][, c("celltype", "donor_id"), drop = FALSE]


# Rebuild object with only counts
pbmc_sc_ref <- CreateSeuratObject(
  counts = counts_mat,
  meta.data = meta_data,
  assay = "RNA",
  project = "pbmc_sc_ref"
)

pbmc_sc_ref@meta.data <- meta_data

# Optional: preserve cell identities if you used celltype as identities
Idents(pbmc_sc_ref) <- pbmc_sc_ref$celltype



pbmc_sc_ref

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  pbmc_sc_ref,
  file = "inst/extdata/pbmc_sc_ref.rds",
  compress = "xz"
)

# Check compressed size
print(file.info("inst/extdata/pbmc_sc_ref.rds")$size / 1024^2)