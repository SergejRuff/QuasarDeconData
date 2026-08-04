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