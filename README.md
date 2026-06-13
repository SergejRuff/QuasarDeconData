
<!-- README.md is generated from README.Rmd. Please edit that file -->

# QuasarDeconData

QuasarDeconData bundles processed peripheral blood mononuclear cell
(PBMC) deconvolution benchmark datasets together with a combined 10x
Genomics single-cell PBMC reference, all formatted for use with QUASAR
deconvolution workflows.

Every dataset shares the same six cell-type labels — `Bcells`,
`CD4Tcells`, `CD8Tcells`, `Monocytes`, `NK`, and `Unknown` — so the
single-cell reference and the bulk benchmarks can be paired directly
without relabeling.

## Installation

You can install the development version of QuasarDeconData from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("SergejRuff/QuasarDeconData")
```

## Datasets

| Object | Platform | Features | Samples / Cells | Ground truth | ID element |
|----|----|----|----|----|----|
| `load_pbmc_sc_ref()` | 10x Genomics scRNA (6k/8k/10k) | 22,629 | 21,978 cells | `celltype` | — |
| `load_cov_imm("train")` | COVID-19 PBMC scRNA | genes in Seurat object | training cells | `cell_type` | `donor_id` |
| `load_cov_imm("test")` | COVID-19 PBMC scRNA | genes in Seurat object | test cells | `cell_type` | `donor_id` |
| `load_cov_pbulk()` | COVID-19 PBMC pseudo-bulk | matched expression features | 10 × 1000 samples | 1000 × 6 per file | — |
| `GSE107011` | RNA-seq | 36,128 | 12 | 12 × 6 | `ensembl_ids` |
| `GSE107572` | RNA-seq | 19,423 | 9 | 9 × 6 | — |
| `GSE120502` | RNA-seq | 39,376 | 249 | 249 × 6 | `gene_ids` |
| `GSE65133` | Microarray (Illumina HT-12 V4) | 35,688 | 20 | 20 × 6 | `probe_ids` |

Note: `GSE65133` is microarray intensity data, not counts — keep that in
mind when feeding it to count-based deconvolution methods.

## Single-cell reference

The reference is a processed Seurat object combining the 10x 6k, 8k, and
10k healthy-donor PBMC datasets. Cell-type labels live in `celltype` and
the dataset of origin in `donor_id`.

``` r
library(QuasarDeconData)

pbmc_sc_ref <- load_pbmc_sc_ref()
pbmc_sc_ref
#> An object of class Seurat 
#> 22629 features across 21978 samples within 1 assay 
#> Active assay: RNA (22629 features, 0 variable features)
#>  1 layer present: counts

# cell-type composition
table(pbmc_sc_ref$celltype)
#> 
#>    Bcells CD4Tcells CD8Tcells Monocytes        NK   Unknown 
#>      3352      8181      3329      5428      1586       102

# dataset of origin
table(pbmc_sc_ref$donor_id)
#> 
#> pbmc10k  pbmc6k  pbmc8k 
#>    9141    4658    8179
```

## Bulk benchmark datasets

Each bulk dataset is a named list holding `bulk_expression_profiles`
(features × samples), `ground_truth_proportions` (samples × cell types),
and a gene/probe identifier vector.

``` r
data(GSE107011)

names(GSE107011)
#> [1] "ground_truth_proportions" "bulk_expression_profiles"
#> [3] "ensembl_ids"

dim(GSE107011$bulk_expression_profiles)
#> [1] 36128    12
dim(GSE107011$ground_truth_proportions)
#> [1] 12  6

# ground-truth proportions sum to 1 per sample
head(GSE107011$ground_truth_proportions)
#>       Bcells CD4Tcells CD8Tcells Monocytes     NK Unknown
#> CYFZ 0.08190   0.17362   0.27818    0.1679 0.1400 0.15840
#> FY2H 0.06927   0.19182   0.17287    0.3267 0.0260 0.21334
#> FLWA 0.08220   0.23634   0.08345    0.2086 0.1950 0.19441
#> 453W 0.05123   0.25331   0.28258    0.2166 0.0678 0.12848
#> 684C 0.09230   0.29714   0.15490    0.2219 0.0845 0.14926
#> CZJE 0.12690   0.19169   0.14691    0.2279 0.1580 0.14860
rowSums(GSE107011$ground_truth_proportions)
#> CYFZ FY2H FLWA 453W 684C CZJE 925L 9JD4 G4YW 4DUY 36TS CR3L 
#>    1    1    1    1    1    1    1    1    1    1    1    1
```

The remaining datasets follow the same structure:

``` r
data(GSE107572)
data(GSE120502)
data(GSE65133)

# samples x cell types for each benchmark
sapply(
  list(
    GSE107011 = GSE107011,
    GSE107572 = GSE107572,
    GSE120502 = GSE120502,
    GSE65133  = GSE65133
  ),
  function(x) dim(x$ground_truth_proportions)
)
#>      GSE107011 GSE107572 GSE120502 GSE65133
#> [1,]        12         9       249       20
#> [2,]         6         6         6        6
```

All four expose the same cell-type columns (column order varies by
dataset):

``` r
lapply(
  list(
    GSE107011 = GSE107011,
    GSE107572 = GSE107572,
    GSE120502 = GSE120502,
    GSE65133  = GSE65133
  ),
  function(x) colnames(x$ground_truth_proportions)
)
#> $GSE107011
#> [1] "Bcells"    "CD4Tcells" "CD8Tcells" "Monocytes" "NK"        "Unknown"  
#> 
#> $GSE107572
#> [1] "CD4Tcells" "CD8Tcells" "NK"        "Bcells"    "Monocytes" "Unknown"  
#> 
#> $GSE120502
#> [1] "Monocytes" "Bcells"    "CD4Tcells" "NK"        "CD8Tcells" "Unknown"  
#> 
#> $GSE65133
#> [1] "Monocytes" "Unknown"   "CD4Tcells" "Bcells"    "NK"        "CD8Tcells"
```

## COVID-19 single-cell train/test reference

QuasarDeconData also includes a donor-balanced COVID-19 PBMC immune
reference split into training and test sets. The split balances mild and
severe COVID-19 donors so that neither half is enriched for disease
severity.

The training split is intended as the single-cell reference for
deconvolution methods, while the test split was used to generate the
packaged COVID-19 pseudo-bulk benchmark datasets.

``` r
train_cov <- load_cov_imm("train")
test_cov  <- load_cov_imm("test")

train_cov
test_cov

# cell-type composition
table(train_cov$cell_type)
table(test_cov$cell_type)

# donor IDs retained for subject-aware methods such as MuSiC
table(train_cov$donor_id)
table(test_cov$donor_id)
```

Both objects are returned as Seurat objects with cell-type identities
set from `cell_type`. The `donor_id` metadata column is retained for
deconvolution methods that require per-subject information.

## COVID-19 pseudo-bulk benchmarks

The COVID-19 pseudo-bulk datasets were generated from `test_cov` and are
intended for deconvolution benchmarking with `train_cov` as the
reference. The package contains ten independent pseudo-bulk files, each
with 1000 simulated bulk samples.

Each pseudo-bulk object contains two elements:

- `bulk_expression_profiles`: pseudo-bulk expression profiles.
- `ground_truth_proportions`: true cell-type proportions for each
  pseudo-bulk sample.

``` r
# load all ten pseudo-bulk datasets
cov_pbulks <- load_cov_pbulk()

length(cov_pbulks)
names(cov_pbulks)

# load one specific pseudo-bulk dataset
cov_pbulk_1 <- load_cov_pbulk(1)

names(cov_pbulk_1)

dim(cov_pbulk_1$bulk_expression_profiles)
dim(cov_pbulk_1$ground_truth_proportions)

# ground-truth proportions sum to 1 per sample
head(cov_pbulk_1$ground_truth_proportions)
rowSums(cov_pbulk_1$ground_truth_proportions)[1:10]
```

A typical benchmark workflow is therefore:

``` r
train_cov <- load_cov_imm("train")
cov_pbulk_1 <- load_cov_pbulk(1)

bulk_expr <- cov_pbulk_1$bulk_expression_profiles
truth <- cov_pbulk_1$ground_truth_proportions

# res = run_deconvolution(bulk = bulk_expr,ref = train_cov)
# metrics = calculate_rsme_pearson(estimates= res,ground_ruth = truth)
```
