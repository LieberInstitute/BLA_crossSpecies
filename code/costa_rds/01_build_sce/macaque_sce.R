library("SingleCellExperiment")
library("DropletUtils")
library("here")
library("rtracklayer")
library("sessioninfo")
library("dplyr")


# directories
raw_dir = here("raw-data")
samples_dir = here(raw_dir, "samples")
sampleinfo_dir = here(raw_dir, "sampleinfo")
processed_dir = here("processed-data")

# ========== Read in the CSV file as a data frame ==========
tmp <- read.delim(here(sampleinfo_dir,"NHP_sampleinfo.csv"), header = T,sep=',')

macaque_info <- tmp %>%
  filter(Species == "Macaque")

# View the subsetted data
head(macaque_info)
# Sample          Region Anatomical Subject Sex Species Batch Seurat.Sample..
# 1     BA1           Basal    Ventral             Macaque     1               1
# 2     BA2           Basal    Ventral             Macaque     1               2
# 3     LA1         Lateral    Ventral             Macaque     1               3
# 4     LA2         Lateral    Ventral             Macaque     1               4
# 5 AB_1A_1 Central Nucleus                        Macaque     1               5
# 6 AB_1A_2 Central Nucleus                        Macaque     1               6

# ========== Set up sample data frame ==========
sample_info <- data.frame(
  #sample_id = paste(tmp$Subject, tmp$Sample, sep = "-"),
  Sample = macaque_info$Sample,
  Region = macaque_info$Region,
  Sex = macaque_info$Sex,
  Anatomoical = macaque_info$Anatomical,
  Batch = macaque_info$Batch
)

head(sample_info)
# Sample          Region Sex Anatomoical Batch
# 1     BA1           Basal         Ventral     1
# 2     BA2           Basal         Ventral     1
# 3     LA1         Lateral         Ventral     1
# 4     LA2         Lateral         Ventral     1
# 5 AB_1A_1 Central Nucleus                     1
# 6 AB_1A_2 Central Nucleus                     1

# check for duplicates
stopifnot(all(!duplicated(sample_info$Sample)))

# add path to cellranger output
sample_info$sample_path
sample_info$sample_path<- file.path(
  samples_dir,
  sample_info$Sample
)

head(sample_info$sample_path)
# [1] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/BA1"    
# [2] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/BA2"    
# [3] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/LA1"    
# [4] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/LA2"    
# [5] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/AB_1A_1"
# [6] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/AB_1A_2"

# ============= Build basic SCE ==============
# I initially tried to run this with both Macacue and Baboon samples and received
# the error "gene information differs between runs".

# After running seperately, it runs fine.
message("Read 10x data and create sce - ", Sys.time())
# Read 10x data and create sce - 2023-08-03 07:42:28

sce.ens <- read10xCounts(
  sample_info$sample_path,
  sample.names = sample_info$Sample,
  type = "sparse",
  col.names = TRUE
)
message("RDone - ", Sys.time())
# RDone - 2023-08-03 07:50:18

# ============= Unique-ify gene names ============= 
# Note that the rownames are the Ensembl gene names - let's switch to the more familiar gene symbols:
# ...but also notice that length(unique(rowData(sce)$Symbol)) != nrow(sce)
#   - That's because some gene symbols are used multiple times.

## Use files from https://ftp.ensembl.org/pub/release-110/gtf/macaca_mulatta/

# gtf <-
#   rtracklayer::import(
#     "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
#   )

gtf <-
   rtracklayer::import(
     "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/genes/Mmul_10.110.gtf"
   )

gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

# Assume 'gtf' is your data frame and 'gene_name' is the column of interest

mt_rows <- subset(gtf, grepl("MT", gene_name))


## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_biotype")]

## Add the gene info to our SPE object
rowRanges(sce) <- gtf[match_genes]

# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq


# ========== Add metadata ===========
# backup.sce <- sce

sce$sample_id <- sce$Sample

sce$key <- paste0(sce$Barcode, "_", sce$Sample)
new_col <- merge(colData(sce), sample_info)
new_col <- new_col[match(sce$key, new_col$key), ]
stopifnot(identical(sce$key, new_col$key))
rownames(new_col) <- colnames(sce)
colData(sce) <- new_col

## Inspect object
sce
# class: SingleCellExperiment 
# dim: 22353 152472 
# metadata(1): Samples
# assays(1): counts
# rownames(22353): U6_ENSMMUG00000036181 ZNF692 ... UBE2M KIR3DH
# rowData names(7): source type ... gene_biotype Symbol.uniq
# colnames(152472): 1_AAACCCAAGAACGCGT-1 1_AAACCCAAGCACCAGA-1 ... 41_TTTGGTTGTGGCTAGA-1 41_TTTGTTGCACATGGTT-1
# colData names(8): Sample Barcode ... Batch sample_path
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):

unique(sce$Batch)
# [1] 1 2 3 4 5

unique(sce$Sample)
# [1] "BA_1"                 "BA_2"                 "LA_1"                 "LA_2"                 "AB_1A_1"             
# [6] "AB_1A_2"              "AB_1B_1"              "AB_1B_2"              "LateralVentralAP1"    "LateralVentralAP2"   
# [11] "BasalDorsalAP1_3A"    "BasalDorsalAP1_3B"    "BasalVentralAP1"      "BasalVentralAP2"      "AccBasalAP1AP2"      
# [16] "LateralDorsalAP1AP2"  "LV_13"                "LV3_13"               "LD_13"                "BV2_13"              
# [21] "BV1_13"               "BD_13"                "AB_13"                "Bd_Bv_13"             "LV2_14"              
# [26] "LV1_14"               "LD_14"                "L_Comb_14"            "BV_14"                "BD_14"               
# [31] "B_Comb_14"            "AB_14"                "CN_S7"                "CN_S8"                "CN_S1"               
# [36] "AccessoryBasalAP2AP3" "BasalAP1_S2"          "BasalAP2_S5"          "Lateral_AP1_S4"       "Lateral_AP2_S6"      
# [41] "Central_Nucleus_AP3" 

unique(sce$Anatomoical)
# [1] "Ventral"        ""               "Dorsal"         "Dorsal_Ventral"

## Size in Gb
lobstr::obj_size(sce) / 1024^3
# 4.44 B

if (!dir.exists(here("processed-data", "01_build_sce"))) dir.create(here("processed-data", "01_build_sce"))
save(sce, file = here("processed-data", "01_build_sce", "sce_maqacue.rda"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.2 (2022-10-31)
# os       macOS Monterey 12.3
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/New_York
# date     2023-08-03
# rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
# pandoc   NA
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version   date (UTC) lib source
# abind                    1.4-5     2016-07-21 [1] CRAN (R 4.2.0)
# AnnotationDbi            1.60.2    2023-03-12 [1] Bioconductor
# AnnotationHub          * 3.6.0     2022-11-07 [1] Bioconductor
# attempt                  0.3.1     2020-05-03 [1] CRAN (R 4.2.0)
# beachmat                 2.14.2    2023-04-09 [1] Bioconductor
# beeswarm                 0.4.0     2021-06-01 [1] CRAN (R 4.2.0)
# benchmarkme              1.0.8     2022-06-12 [1] CRAN (R 4.2.0)
# benchmarkmeData          1.0.4     2020-04-23 [1] CRAN (R 4.2.0)
# Biobase                * 2.58.0    2022-11-07 [1] Bioconductor
# BiocFileCache          * 2.6.1     2023-02-19 [1] Bioconductor
# BiocGenerics           * 0.44.0    2022-11-07 [1] Bioconductor
# BiocIO                   1.8.0     2022-11-07 [1] Bioconductor
# BiocManager              1.30.20   2023-02-24 [1] CRAN (R 4.2.0)
# BiocNeighbors            1.16.0    2022-11-07 [1] Bioconductor
# BiocParallel             1.32.6    2023-03-19 [1] Bioconductor
# BiocSingular             1.14.0    2022-11-07 [1] Bioconductor
# BiocVersion              3.16.0    2022-09-20 [1] Bioconductor
# Biostrings               2.66.0    2022-11-07 [1] Bioconductor
# bit                      4.0.5     2022-11-15 [1] CRAN (R 4.2.0)
# bit64                    4.0.5     2020-08-30 [1] CRAN (R 4.2.0)
# bitops                   1.0-7     2021-04-24 [1] CRAN (R 4.2.0)
# blob                     1.2.4     2023-03-17 [1] CRAN (R 4.2.0)
# bluster                  1.8.0     2022-11-07 [1] Bioconductor
# bslib                    0.4.2     2022-12-16 [1] CRAN (R 4.2.0)
# cachem                   1.0.8     2023-05-01 [1] CRAN (R 4.2.0)
# cli                      3.6.1     2023-03-23 [1] CRAN (R 4.2.0)
# cluster                  2.1.4     2022-08-22 [1] CRAN (R 4.2.2)
# codetools                0.2-19    2023-02-01 [1] CRAN (R 4.2.0)
# colorspace               2.1-0     2023-01-23 [1] CRAN (R 4.2.0)
# config                   0.3.1     2020-12-17 [1] CRAN (R 4.2.0)
# cowplot                  1.1.1     2020-12-30 [1] CRAN (R 4.2.0)
# crayon                   1.5.2     2022-09-29 [1] CRAN (R 4.2.0)
# curl                     5.0.0     2023-01-12 [1] CRAN (R 4.2.0)
# data.table               1.14.8    2023-02-17 [1] CRAN (R 4.2.0)
# DBI                      1.1.3     2022-06-18 [1] CRAN (R 4.2.0)
# dbplyr                 * 2.3.2     2023-03-21 [1] CRAN (R 4.2.0)
# DelayedArray             0.24.0    2022-11-07 [1] Bioconductor
# DelayedMatrixStats       1.20.0    2022-11-01 [1] Bioconductor
# deldir                   1.0-6     2021-10-23 [1] CRAN (R 4.2.0)
# digest                   0.6.31    2022-12-11 [1] CRAN (R 4.2.0)
# doParallel               1.0.17    2022-02-07 [1] CRAN (R 4.2.0)
# dotCall64                1.0-2     2022-10-03 [1] CRAN (R 4.2.0)
# dplyr                  * 1.1.2     2023-04-20 [1] CRAN (R 4.2.0)
# dqrng                    0.3.0     2021-05-01 [1] CRAN (R 4.2.0)
# DropletUtils           * 1.18.1    2022-11-23 [1] Bioconductor
# DT                       0.27      2023-01-17 [1] CRAN (R 4.2.0)
# edgeR                    3.40.2    2023-01-22 [1] Bioconductor
# ellipsis                 0.3.2     2021-04-29 [1] CRAN (R 4.2.0)
# ExperimentHub            2.6.0     2022-11-07 [1] Bioconductor
# fansi                    1.0.4     2023-01-22 [1] CRAN (R 4.2.0)
# fastmap                  1.1.1     2023-02-24 [1] CRAN (R 4.2.0)
# fields                   14.1      2022-08-12 [1] CRAN (R 4.2.0)
# filelock                 1.0.2     2018-10-05 [1] CRAN (R 4.2.0)
# fitdistrplus             1.1-11    2023-04-25 [1] CRAN (R 4.2.0)
# foreach                  1.5.2     2022-02-02 [1] CRAN (R 4.2.0)
# future                   1.32.0    2023-03-07 [1] CRAN (R 4.2.0)
# future.apply             1.10.0    2022-11-05 [1] CRAN (R 4.2.0)
# generics                 0.1.3     2022-07-05 [1] CRAN (R 4.2.0)
# GenomeInfoDb           * 1.34.9    2023-02-05 [1] Bioconductor
# GenomeInfoDbData         1.2.9     2022-11-29 [1] Bioconductor
# GenomicAlignments        1.34.1    2023-03-12 [1] Bioconductor
# GenomicRanges          * 1.50.2    2022-12-18 [1] Bioconductor
# ggbeeswarm               0.7.2     2023-04-29 [1] CRAN (R 4.2.0)
# ggplot2                * 3.4.2     2023-04-03 [1] CRAN (R 4.2.0)
# ggrepel                  0.9.3     2023-02-03 [1] CRAN (R 4.2.0)
# ggridges                 0.5.4     2022-09-26 [1] CRAN (R 4.2.0)
# ggside                   0.2.2     2022-12-04 [1] CRAN (R 4.2.0)
# ggspavis               * 1.4.0     2022-11-07 [1] Bioconductor
# globals                  0.16.2    2022-11-21 [1] CRAN (R 4.2.0)
# glue                     1.6.2     2022-02-24 [1] CRAN (R 4.2.0)
# goftest                  1.2-3     2021-10-07 [1] CRAN (R 4.2.0)
# golem                    0.4.0     2023-03-12 [1] CRAN (R 4.2.0)
# gridExtra                2.3       2017-09-09 [1] CRAN (R 4.2.0)
# gtable                   0.3.3     2023-03-21 [1] CRAN (R 4.2.0)
# HDF5Array                1.26.0    2022-11-07 [1] Bioconductor
# here                   * 1.0.1     2020-12-13 [1] CRAN (R 4.2.0)
# htmltools                0.5.5     2023-03-23 [1] CRAN (R 4.2.0)
# htmlwidgets              1.6.2     2023-03-17 [1] CRAN (R 4.2.0)
# httpuv                   1.6.11    2023-05-11 [1] CRAN (R 4.2.2)
# httr                     1.4.6     2023-05-08 [1] CRAN (R 4.2.2)
# ica                      1.0-3     2022-07-08 [1] CRAN (R 4.2.0)
# igraph                   1.4.2     2023-04-07 [1] CRAN (R 4.2.0)
# interactiveDisplayBase   1.36.0    2022-11-07 [1] Bioconductor
# IRanges                * 2.32.0    2022-11-07 [1] Bioconductor
# irlba                    2.3.5.1   2022-10-03 [1] CRAN (R 4.2.0)
# iterators                1.0.14    2022-02-05 [1] CRAN (R 4.2.0)
# jquerylib                0.1.4     2021-04-26 [1] CRAN (R 4.2.0)
# jsonlite                 1.8.4     2022-12-06 [1] CRAN (R 4.2.2)
# KEGGREST                 1.38.0    2022-11-07 [1] Bioconductor
# KernSmooth               2.23-21   2023-05-03 [1] CRAN (R 4.2.0)
# later                    1.3.1     2023-05-02 [1] CRAN (R 4.2.0)
# lattice                  0.21-8    2023-04-05 [1] CRAN (R 4.2.0)
# lazyeval                 0.2.2     2019-03-15 [1] CRAN (R 4.2.0)
# leiden                   0.4.3     2022-09-10 [1] CRAN (R 4.2.0)
# lifecycle                1.0.3     2022-10-07 [1] CRAN (R 4.2.0)
# limma                    3.54.2    2023-03-01 [1] Bioconductor
# listenv                  0.9.0     2022-12-16 [1] CRAN (R 4.2.0)
# lmtest                   0.9-40    2022-03-21 [1] CRAN (R 4.2.0)
# lobstr                   1.1.2     2022-06-22 [1] CRAN (R 4.2.0)
# locfit                   1.5-9.7   2023-01-02 [1] CRAN (R 4.2.2)
# magick                   2.7.4     2023-03-09 [1] CRAN (R 4.2.0)
# magrittr                 2.0.3     2022-03-30 [1] CRAN (R 4.2.0)
# maps                     3.4.1     2022-10-30 [1] CRAN (R 4.2.0)
# MASS                     7.3-60    2023-05-04 [1] CRAN (R 4.2.0)
# Matrix                 * 1.5-4     2023-04-04 [1] CRAN (R 4.2.0)
# MatrixGenerics         * 1.10.0    2022-11-07 [1] Bioconductor
# matrixStats            * 0.63.0    2022-11-18 [1] CRAN (R 4.2.0)
# memoise                  2.0.1     2021-11-26 [1] CRAN (R 4.2.0)
# metapod                  1.6.0     2022-11-07 [1] Bioconductor
# mime                     0.12      2021-09-28 [1] CRAN (R 4.2.0)
# miniUI                   0.1.1.1   2018-05-18 [1] CRAN (R 4.2.0)
# munsell                  0.5.0     2018-06-12 [1] CRAN (R 4.2.0)
# nlme                     3.1-162   2023-01-31 [1] CRAN (R 4.2.0)
# paletteer                1.5.0     2022-10-19 [1] CRAN (R 4.2.0)
# parallelly               1.35.0    2023-03-23 [1] CRAN (R 4.2.0)
# patchwork              * 1.1.2     2022-08-19 [1] CRAN (R 4.2.0)
# pbapply                  1.7-0     2023-01-13 [1] CRAN (R 4.2.0)
# pheatmap               * 1.0.12    2019-01-04 [1] CRAN (R 4.2.0)
# pillar                   1.9.0     2023-03-22 [1] CRAN (R 4.2.0)
# pkgconfig                2.0.3     2019-09-22 [1] CRAN (R 4.2.0)
# plotly                   4.10.1    2022-11-07 [1] CRAN (R 4.2.0)
# plyr                     1.8.8     2022-11-11 [1] CRAN (R 4.2.0)
# png                      0.1-8     2022-11-29 [1] CRAN (R 4.2.0)
# polyclip                 1.10-4    2022-10-20 [1] CRAN (R 4.2.0)
# prettyunits              1.1.1     2020-01-24 [1] CRAN (R 4.2.0)
# progressr                0.13.0    2023-01-10 [1] CRAN (R 4.2.0)
# promises                 1.2.0.1   2021-02-11 [1] CRAN (R 4.2.0)
# purrr                    1.0.1     2023-01-10 [1] CRAN (R 4.2.0)
# R.methodsS3              1.8.2     2022-06-13 [1] CRAN (R 4.2.0)
# R.oo                     1.25.0    2022-06-12 [1] CRAN (R 4.2.0)
# R.utils                  2.12.2    2022-11-11 [1] CRAN (R 4.2.0)
# R6                       2.5.1     2021-08-19 [1] CRAN (R 4.2.0)
# RANN                     2.6.1     2019-01-08 [1] CRAN (R 4.2.0)
# rappdirs                 0.3.3     2021-01-31 [1] CRAN (R 4.2.0)
# RColorBrewer             1.1-3     2022-04-03 [1] CRAN (R 4.2.0)
# Rcpp                     1.0.10    2023-01-22 [1] CRAN (R 4.2.0)
# RcppAnnoy                0.0.20    2022-10-27 [1] CRAN (R 4.2.0)
# RCurl                    1.98-1.12 2023-03-27 [1] CRAN (R 4.2.0)
# rematch2                 2.1.2     2020-05-01 [1] CRAN (R 4.2.0)
# reshape2                 1.4.4     2020-04-09 [1] CRAN (R 4.2.0)
# restfulr                 0.0.15    2022-06-16 [1] CRAN (R 4.2.0)
# reticulate               1.28      2023-01-27 [1] CRAN (R 4.2.0)
# rhdf5                    2.42.1    2023-04-09 [1] Bioconductor
# rhdf5filters             1.10.1    2023-03-26 [1] Bioconductor
# Rhdf5lib                 1.20.0    2022-11-07 [1] Bioconductor
# rjson                    0.2.21    2022-01-09 [1] CRAN (R 4.2.0)
# rlang                    1.1.1     2023-04-28 [1] CRAN (R 4.2.0)
# ROCR                     1.0-11    2020-05-02 [1] CRAN (R 4.2.0)
# rprojroot                2.0.3     2022-04-02 [1] CRAN (R 4.2.0)
# Rsamtools                2.14.0    2022-11-07 [1] Bioconductor
# RSQLite                  2.3.1     2023-04-03 [1] CRAN (R 4.2.0)
# rstudioapi               0.14      2022-08-22 [1] CRAN (R 4.2.0)
# rsvd                     1.0.5     2021-04-16 [1] CRAN (R 4.2.0)
# rtracklayer            * 1.58.0    2022-11-07 [1] Bioconductor
# Rtsne                    0.16      2022-04-17 [1] CRAN (R 4.2.0)
# S4Vectors              * 0.36.2    2023-03-01 [1] Bioconductor
# sass                     0.4.6     2023-05-03 [1] CRAN (R 4.2.0)
# ScaledMatrix             1.6.0     2022-11-07 [1] Bioconductor
# scales                   1.2.1     2022-08-20 [1] CRAN (R 4.2.0)
# scater                 * 1.26.1    2022-11-13 [1] Bioconductor
# scattermore              1.0       2023-05-03 [1] CRAN (R 4.2.0)
# scran                  * 1.26.2    2023-01-22 [1] Bioconductor
# sctransform              0.3.5     2022-09-21 [1] CRAN (R 4.2.0)
# scuttle                * 1.8.4     2023-01-22 [1] Bioconductor
# sessioninfo            * 1.2.2     2021-12-06 [1] CRAN (R 4.2.0)
# Seurat                 * 4.3.0     2022-11-18 [1] CRAN (R 4.2.0)
# SeuratObject           * 4.1.3     2022-11-07 [1] CRAN (R 4.2.0)
# shiny                    1.7.4     2022-12-15 [1] CRAN (R 4.2.0)
# shinyWidgets             0.7.6     2023-01-08 [1] CRAN (R 4.2.0)
# SingleCellExperiment   * 1.20.1    2023-03-19 [1] Bioconductor
# sp                       1.6-0     2023-01-19 [1] CRAN (R 4.2.0)
# spam                     2.9-1     2022-08-07 [1] CRAN (R 4.2.0)
# sparseMatrixStats        1.10.0    2022-11-07 [1] Bioconductor
# SpatialExperiment      * 1.8.1     2023-03-05 [1] Bioconductor
# spatialLIBD            * 1.13.3    2023-05-12 [1] Github (LieberInstitute/spatialLIBD@7a3d71d)
# spatstat.data            3.0-1     2023-03-12 [1] CRAN (R 4.2.0)
# spatstat.explore         3.1-0     2023-03-14 [1] CRAN (R 4.2.0)
# spatstat.geom            3.2-1     2023-05-09 [1] CRAN (R 4.2.2)
# spatstat.random          3.1-5     2023-05-11 [1] CRAN (R 4.2.2)
# spatstat.sparse          3.0-1     2023-03-12 [1] CRAN (R 4.2.0)
# spatstat.utils           3.0-3     2023-05-09 [1] CRAN (R 4.2.2)
# statmod                  1.5.0     2023-01-06 [1] CRAN (R 4.2.0)
# stringi                  1.7.12    2023-01-11 [1] CRAN (R 4.2.0)
# stringr                  1.5.0     2022-12-02 [1] CRAN (R 4.2.0)
# SummarizedExperiment   * 1.28.0    2022-11-07 [1] Bioconductor
# survival                 3.5-5     2023-03-12 [1] CRAN (R 4.2.0)
# tensor                   1.5       2012-05-05 [1] CRAN (R 4.2.0)
# tibble                   3.2.1     2023-03-20 [1] CRAN (R 4.2.0)
# tidyr                    1.3.0     2023-01-24 [1] CRAN (R 4.2.0)
# tidyselect               1.2.0     2022-10-10 [1] CRAN (R 4.2.0)
# utf8                     1.2.3     2023-01-31 [1] CRAN (R 4.2.0)
# uwot                     0.1.14    2022-08-22 [1] CRAN (R 4.2.0)
# vctrs                    0.6.3     2023-06-14 [1] CRAN (R 4.2.0)
# vipor                    0.4.5     2017-03-22 [1] CRAN (R 4.2.0)
# viridis                  0.6.3     2023-05-03 [1] CRAN (R 4.2.0)
# viridisLite              0.4.2     2023-05-02 [1] CRAN (R 4.2.0)
# withr                    2.5.0     2022-03-03 [1] CRAN (R 4.2.0)
# XML                      3.99-0.14 2023-03-19 [1] CRAN (R 4.2.0)
# xtable                   1.8-4     2019-04-21 [1] CRAN (R 4.2.0)
# XVector                  0.38.0    2022-11-07 [1] Bioconductor
# yaml                     2.3.7     2023-01-23 [1] CRAN (R 4.2.0)
# zlibbioc                 1.44.0    2022-11-07 [1] Bioconductor
# zoo                      1.8-12    2023-04-13 [1] CRAN (R 4.2.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# > 

  