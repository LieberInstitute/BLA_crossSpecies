# qrsh -l mem_free=80G,h_vmem=80G

library("SingleCellExperiment")
library("DropletUtils")
library("here")
library("rtracklayer")
library("sessioninfo")
library("dplyr")

# Read in the CSV file as a data frame
tmp <- read.delim(here("raw-data","sampleinfo",
                "master_sampleinfo_2023-10-09.csv"),
                header = T,sep=',')

# View the subsetted data
head(tmp)

#        Sample Species Subject Sex   Region Subregion DV_axis  PI.NeuN
#   1 3c-AMYBLA   Human  Br2327     Amygdala       BLA         PI+NeuN+
#   2 4c-AMYBLA   Human  Br8692     Amygdala       BLA         PI+NeuN+
#   3 5c-AMYBLA   Human  Br9021     Amygdala       BLA         PI+NeuN+
#   4  34ac_scp   Human  Br8331     Amygdala       BLA         PI+NeuN+
#   5  35ac_scp   Human  Br5273     Amygdala       BLA         PI+NeuN+
#   6   Sample1 Macaque    Mac2     Amygdala   Lateral Ventral         

##set up sample data table
sample_info <- data.frame(
  sample_id = paste(tmp$Subject, tmp$Sample, sep = "-"),
  sample_name = tmp$Sample,
  subject = tmp$Subject,
  species = tmp$Species,
  region = tmp$Region,
  subregion = tmp$Subregion,
  dv_axis = tmp$DV_axis
)

stopifnot(all(!duplicated(sample_info$sample_id)))

# add path to cellranger output
sample_info$sample_path<-rep(NA,2)
sample_info$sample_path<- file.path(
  here::here("processed-data", "01_cellranger"),
  sample_info$species,
  sample_info$sample_name,
  "outs",
  "raw_feature_bc_matrix"
)

# Subset to just human samples
sample_info <- sample_info %>%
  filter(species == "Human")


## Build basic SCE
message("Read 10x data and create sce - ", Sys.time())
# Read 10x data and create sce - 2023-04-17 13:37:16

sce <- read10xCounts(
  sample_info$sample_path,
  sample_info$sample_id,
  type = "sparse",
  col.names = TRUE
)
message("RDone - ", Sys.time())
# RDone - 2023-04-17 13:38:35

subset_sce <- sce[, sce$Sample == "5c-AMYBLA"]

# Note that the rownames are the Ensembl gene names - let's switch to the more familiar gene symbols:
# ...but also notice that length(unique(rowData(sce)$Symbol)) != nrow(sce)
#   - That's because some gene symbols are used multiple times.

## Use code from https://github.com/LieberInstitute/Visium_IF_AD/commit/08df3f7e4a3178563d6b4b1861b664b21466b395#diff-10cb35de98e2a3e5f4235cd88f6dabce5469eead2b2db1fd7121126849fcf585L100
## Read in the gene information from the annotation GTF file
gtf <-
    rtracklayer::import(
        "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our SPE object
rowRanges(sce) <- gtf[match_genes]

# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq


# Add metadata
sce$key <- paste0(sce$Barcode, "_", sce$Sample)
new_col <- merge(colData(sce), sample_info[, -which(colnames(sample_info) == "sample_path")])
new_col$key <- paste0(new_col$Barcode, "_", new_col$Sample)
new_col <- new_col[match(sce$key, new_col$key), ]
stopifnot(identical(sce$key, new_col$key))
rownames(new_col) <- colnames(sce)
colData(sce) <- new_col

## Inspect object
sce
# class: SingleCellExperiment 
# dim: 36601 7832046 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(7832046): 1_AAACCCAAGAAACACT-1 1_AAACCCAAGAAACCCA-1 ...
# 5_TTTGTTGTCTTTGATC-1 5_TTTGTTGTCTTTGCTA-1
# colData names(10): Sample Barcode ... dv_axis key
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):


if (!dir.exists(here("processed-data", "02_build_sce"))) dir.create(here("processed-data", "02_build_sce"))
save(sce, file = here("processed-data", "02_build_sce", "sce_human_raw.rda"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.3 Patched (2023-04-07 r84211)
# os       CentOS Linux 7 (Core)
# m
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.2.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.2.3)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.2.3)
# colorout             * 1.2-2     2023-02-15 [1] Github (jalvesaq/colorout@79931fd)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.2.1)
# DelayedArray           0.24.0    2022-11-01 [2] Bioconductor
# DelayedMatrixStats     1.20.0    2022-11-01 [2] Bioconductor
# dplyr                * 1.1.1     2023-03-22 [2] CRAN (R 4.2.3)
# dqrng                  0.3.0     2021-05-01 [2] CRAN (R 4.2.1)
# DropletUtils         * 1.18.1    2022-11-22 [2] Bioconductor
# edgeR                  3.40.2    2023-01-19 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.2.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.2.1)
# GenomeInfoDb         * 1.34.9    2023-02-02 [2] Bioconductor
# GenomeInfoDbData       1.2.9     2022-09-29 [2] Bioconductor
# GenomicAlignments      1.34.1    2023-03-09 [2] Bioconductor
# GenomicRanges        * 1.50.2    2022-12-16 [2] Bioconductor
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.2.1)
# HDF5Array              1.26.0    2022-11-01 [2] Bioconductor
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.2.1)
# IRanges              * 2.32.0    2022-11-01 [2] Bioconductor
# lattice                0.20-45   2021-09-22 [3] CRAN (R 4.2.3)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.2.1)
# limma                  3.54.2    2023-02-28 [2] Bioconductor
# locfit                 1.5-9.7   2023-01-02 [2] CRAN (R 4.2.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.2.1)
# Matrix                 1.5-4     2023-04-04 [2] CRAN (R 4.2.3)
# MatrixGenerics       * 1.10.0    2022-11-01 [2] Bioconductor
# matrixStats          * 0.63.0    2022-11-18 [2] CRAN (R 4.2.2)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.2.3)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.2.1)
# R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.2.1)
# R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.2.1)
# R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.2.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.2.1)
# Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.2.2)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.2.3)
# restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.2.1)
# rhdf5                  2.42.1    2023-04-07 [2] Bioconductor
# rhdf5filters           1.10.1    2023-03-24 [2] Bioconductor
# Rhdf5lib               1.20.0    2022-11-01 [2] Bioconductor
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.2.1)
# rlang                  1.1.0     2023-03-14 [2] CRAN (R 4.2.3)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.2.1)
# Rsamtools              2.14.0    2022-11-01 [2] Bioconductor
# rtracklayer          * 1.58.0    2022-11-01 [2] Bioconductor
# S4Vectors            * 0.36.2    2023-02-26 [2] Bioconductor
# scuttle                1.8.4     2023-01-19 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.2.1)
# SingleCellExperiment * 1.20.1    2023-03-17 [2] Bioconductor
# sparseMatrixStats      1.10.0    2022-11-01 [2] Bioconductor
# SummarizedExperiment * 1.28.0    2022-11-01 [2] Bioconductor
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.2.3)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.2.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.2.2)
# vctrs                  0.6.1     2023-03-22
