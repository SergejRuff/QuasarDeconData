
<!-- README.md is generated from README.Rmd. Please edit that file -->

# QuasarDeconData

QuasarDeconData bundles processed peripheral blood mononuclear cell
(PBMC) deconvolution benchmark datasets together with a combined 10x
Genomics single-cell PBMC reference, all formatted for use with QUASAR
deconvolution workflows.

The bulk benchmark datasets listed below and the healthy 10x PBMC
reference use the same six cell-type labels in their harmonized
annotations: `Bcells`, `CD4Tcells`, `CD8Tcells`, `Monocytes`, `NK`, and
`Unknown`. Column order and gene identifier formats vary, so align them
before benchmarking.

The COVID-19 references and pseudo-bulk datasets retain their own
cell-type annotations.

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
| `load_cov_imm("healthy")` | COVID-19 PBMC scRNA (healthy donors only) | genes in Seurat object | healthy cells | `cell_type` | `donor_id` |
| `load_cov_pbulk()` | COVID-19 PBMC pseudo-bulk | matched expression features | 10 × 1000 samples | 1000 × 11 per file | — |
| `load_cov_pbulk("train")` | COVID-19 PBMC pseudo-bulk | matched expression features | 10,000 samples | 10,000 × 11 | — |
| `load_cov_pbulk("healthy")` | COVID-19 PBMC pseudo-bulk (healthy donors only) | matched expression features | healthy samples | matched proportions | — |
| `GSE107011` (Monaco1) | RNA-seq | 36,128 | 12 | 12 × 6 | `ensembl_ids` |
| `GSE107990` (Monaco2) | Microarray (Illumina HT-12 V4) | 20,908 | 164 | 164 × 6; additional 164 × 11 | — |
| `GSE107572` | RNA-seq | 19,423 | 9 | 9 × 6 | — |
| `GSE120502` | RNA-seq | 39,376 | 249 | 249 × 6 | `gene_ids` |
| `GSE65133` | Microarray (Illumina HT-12 V4) | 35,688 | 20 | 20 × 6 | `probe_ids` |
| `GSE193141` | RNA-seq (raw counts) | 39,376 | 156 | 156 × 6; additional 156 × 8 | — |
| `SDY67` | RNA-seq | 23,113 | 231 | 231 × 6 | — |
| `SDY67_scaden` | RNA-seq | 11,328 | 12 | 12 × 6 | — |

Note: `GSE107990` and `GSE65133` contain microarray intensity
measurements rather than sequencing counts.

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

Each bulk dataset is a named list containing:

- `bulk_expression_profiles`: expression values arranged as features ×
  samples.
- `ground_truth_proportions`: reference proportions arranged as samples
  × cell types.

Some datasets also contain a separate gene or probe identifier vector.
For `SDY67` and `GSE107990`, gene symbols are stored in the expression
row names. For `GSE193141`, the expression row names contain NCBI gene
identifiers.

Two datasets additionally contain `ground_truth_all_ct`: `GSE193141`
retains eight original population measurements, and `GSE107990` retains
eleven.

``` r
data(GSE107011)

names(GSE107011)
#> [1] "ground_truth_proportions" "bulk_expression_profiles"
#> [3] "ensembl_ids"

dim(GSE107011$bulk_expression_profiles)
#> [1] 36128    12
dim(GSE107011$ground_truth_proportions)
#> [1] 12  6

# Six-column reference proportions
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

Load the remaining datasets and compare their dimensions:

``` r
data(GSE107990)
data(GSE107572)
data(GSE120502)
data(GSE65133)
data(SDY67)
data(GSE193141)
data(SDY67_scaden)
data(GSE107011)

bulk_datasets <- list(
  GSE107011    = GSE107011,
  GSE107990    = GSE107990,
  GSE107572    = GSE107572,
  GSE120502    = GSE120502,
  GSE65133     = GSE65133,
  SDY67        = SDY67,
  SDY67_scaden = SDY67_scaden,
  GSE193141    = GSE193141
)

t(vapply(
  bulk_datasets,
  function(x) {
    c(
      features = nrow(x$bulk_expression_profiles),
      samples = ncol(x$bulk_expression_profiles),
      cell_types = ncol(x$ground_truth_proportions)
    )
  },
  integer(3)
))
#>              features samples cell_types
#> GSE107011       36128      12          6
#> GSE107990       20908     164          6
#> GSE107572       19423       9          6
#> GSE120502       39376     249          6
#> GSE65133        21186      20          6
#> SDY67           23113     231          6
#> SDY67_scaden    11328      12          6
#> GSE193141       39376     156          6
```

All listed bulk benchmarks expose the same cell-type labels in
`ground_truth_proportions`, although their column order varies:

``` r
lapply(
  bulk_datasets,
  function(x) colnames(x$ground_truth_proportions)
)
#> $GSE107011
#> [1] "Bcells"    "CD4Tcells" "CD8Tcells" "Monocytes" "NK"        "Unknown"  
#> 
#> $GSE107990
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
#> 
#> $SDY67
#> [1] "CD8Tcells" "Monocytes" "CD4Tcells" "Bcells"    "NK"        "Unknown"  
#> 
#> $SDY67_scaden
#> [1] "Monocytes" "Unknown"   "CD4Tcells" "Bcells"    "NK"        "CD8Tcells"
#> 
#> $GSE193141
#> [1] "CD8Tcells" "Monocytes" "CD4Tcells" "Bcells"    "NK"        "Unknown"
```

### SDY67

`SDY67` contains 23,113 genes and 231 samples with complete six-column
proportions. RNA and flow cytometry measurements were matched by
participant and collection visit.

The proportions combine measurements from the innate and Treg panels:

- `Monocytes` combines classical, intermediate, and non-classical
  monocytes.
- `CD4Tcells`, `Bcells`, and `NK` use their corresponding source gates.
- `CD8Tcells` uses the study-labeled CD8 gate, which is a
  CD4-low/negative T-cell proxy. CD8 was not directly measured by this
  panel.
- `Unknown` is calculated as one minus the sum of the five selected
  populations.

The combined proportions are an approximate composition. `Unknown` can
include unassigned cells and differences between the panels’ gates.
Samples with missing or ambiguous panel matches, or overlapping main
CD4/CD8-proxy gates, were excluded.

`SDY67_scaden` are the same 12 bulk samples along the corresponding
ground truth proportion vectors used in the scaden paper and other works
such the DISSECT paper.

``` r
dim(SDY67$bulk_expression_profiles)
#> [1] 23113   231
dim(SDY67$ground_truth_proportions)
#> [1] 231   6

# Expression columns and proportion rows have identical sample ordering
identical(
  colnames(SDY67$bulk_expression_profiles),
  rownames(SDY67$ground_truth_proportions)
)
#> [1] TRUE
```

## COVID-19 single-cell references

QuasarDeconData also includes a donor-balanced COVID-19 PBMC immune
reference split into training and test sets. The split balances mild and
severe COVID-19 donors so that neither half is enriched for disease
severity.

The training split is intended as the single-cell reference for
deconvolution methods, while the test split was used to generate the
packaged COVID-19 pseudo-bulk benchmark datasets.

A third reference, `load_cov_imm("healthy")`, contains cells from
non-diseased donors only. It exists for reference-mismatch experiments,
where a reference that does not match the disease composition of the
mixtures is deliberately used for deconvolution.

``` r
train_cov   <- load_cov_imm("train")
test_cov    <- load_cov_imm("test")
healthy_cov <- load_cov_imm("healthy")

train_cov
test_cov
healthy_cov

# cell-type composition
table(train_cov$cell_type)
table(test_cov$cell_type)
table(healthy_cov$cell_type)

# donor IDs retained for subject-aware methods such as MuSiC
table(train_cov$donor_id)
table(test_cov$donor_id)
table(healthy_cov$donor_id)
```

All three objects are returned as Seurat objects with cell-type
identities set from `cell_type`. The `donor_id` metadata column is
retained for deconvolution methods that require per-subject information.

Note that the healthy reference may lack cell types that are present in
the pseudo-bulk ground truth. Methods cannot recover a cell type absent
from the reference, so predicted proportion matrices from the healthy
reference can have fewer columns than the corresponding
`ground_truth_proportions`. Align the two explicitly — for example by
padding the missing cell types with zeros — before computing benchmark
metrics.

## COVID-19 pseudo-bulk benchmarks

The COVID-19 pseudo-bulk datasets were generated from `test_cov` and are
intended for deconvolution benchmarking with `train_cov` as the
reference. The package contains ten independent pseudo-bulk files, each
with 1000 simulated bulk samples.

Pseudo-bulk data generated from the training split (10,000 samples) and
from the healthy-only split are also available. Both are stored as
chunks internally and reassembled on load, so the split is not visible
to the caller.

Each pseudo-bulk object contains two elements:

- `bulk_expression_profiles`: pseudo-bulk expression profiles.
- `ground_truth_proportions`: true cell-type proportions for each
  pseudo-bulk sample.

``` r
# load all ten test pseudo-bulk datasets
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

# pseudo-bulk from the training split, and from the healthy-only split
train_pbulk   <- load_cov_pbulk("train")
healthy_pbulk <- load_cov_pbulk("healthy")

dim(train_pbulk$bulk_expression_profiles)
dim(healthy_pbulk$bulk_expression_profiles)
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

### Reference-mismatch benchmarking

Swapping the matched reference for the healthy-only one gives a
benchmark in which the reference no longer reflects the disease
composition of the mixtures. Comparing the two runs on the same
pseudo-bulk data isolates the effect of that mismatch on each method:

``` r
cov_pbulk_1 <- load_cov_pbulk(1)

matched  <- load_cov_imm("train")     # reference matching the mixtures
mismatch <- load_cov_imm("healthy")   # healthy-donor reference

# cell types the healthy reference cannot recover
setdiff(
  colnames(cov_pbulk_1$ground_truth_proportions),
  unique(as.character(mismatch$cell_type))
)

# res_matched  = run_deconvolution(bulk = cov_pbulk_1$bulk_expression_profiles, ref = matched)
# res_mismatch = run_deconvolution(bulk = cov_pbulk_1$bulk_expression_profiles, ref = mismatch)
```
