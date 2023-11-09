library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)

## save directories
plot_dir = here("plots", "04_normalization")
processed_dir = here("processed-data","04_normalization")

# load sce
load(file=here("processed-data","03_quality_control","Human","sce_post_qc.rda"))
sce.human <- sce
sce.human
# class: SingleCellExperiment 
# dim: 36601 21268 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(21268): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 5_TTTGTTGTCGGACTGC-1 5_TTTGTTGTCGTTGTTT-1
# colData names(27): Sample Barcode ... discard_auto doubletScore
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

load(file=here("processed-data","03_quality_control","Macaque","PerCellQC","sce_post_qc.rda"))
sce.macaque <- sce
sce.macaque
# class: SingleCellExperiment 
# dim: 21369 109142 
# metadata(1): Samples
# assays(1): counts
# rownames(21369): ENSMMUG00000023296 ZNF692 ... ND6 CYTB
# rowData names(7): source type ... gene_biotype Symbol.uniq
# colnames(109142): 1_AAACGAAAGCGCCGTT-1 1_AAACGAACAGCCGTTG-1 ...
# 35_TTTGGTTGTGGCTAGA-1 35_TTTGTTGCACATGGTT-1
# colData names(13): Sample Barcode ... discard_auto doubletScore
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

load(file=here("processed-data","03_quality_control","Baboon","PerCellQC","sce_post_qc.rda"))
sce.baboon <- sce
sce.baboon
# class: SingleCellExperiment 
# dim: 20842 57564 
# metadata(1): Samples
# assays(1): counts
# rownames(20842): PITHD1 ACTL8 ... ENSPANG00000037790 ENSPANG00000045064
# rowData names(7): source type ... gene_biotype Symbol.uniq
# colnames(57564): 1_AAACCCAAGACTGTTC-1 1_AAACCCAAGATCCAAA-1 ...
# 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# colData names(21): Sample Barcode ... discard_auto doubletScore
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

# =========== Calculate size factors ==============

lib.sf.human <- librarySizeFactors(sce.human)
lib.sf.macaque <- librarySizeFactors(sce.macaque)
lib.sf.baboon <- librarySizeFactors(sce.baboon)

summary(lib.sf.human)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01228 0.41602 0.85299 1.00000 1.40823 5.98530 

summary(lib.sf.macaque)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.02374  0.22615  0.53045  1.00000  1.35222 35.98417

summary(lib.sf.baboon)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.07627  0.40058  0.75330  1.00000  1.29658 14.97605 


# plot size factos
pdf(here(plot_dir,"size_factors.pdf"), width=8, height=4)
hist(log10(lib.sf.human), xlab="Log10[Size factor]", col='grey80')
hist(log10(lib.sf.macaque), xlab="Log10[Size factor]", col='grey80')
hist(log10(lib.sf.baboon), xlab="Log10[Size factor]", col='grey80')
dev.off()

# ========= Normaliz based on size factors ========

sce.human <- logNormCounts(sce.human)
sce.human
# class: SingleCellExperiment 
# dim: 36601 21268 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(21268): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 5_TTTGTTGTCGGACTGC-1 5_TTTGTTGTCGTTGTTT-1
# colData names(28): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

sce.macaque <- logNormCounts(sce.macaque)
sce.macaque
# class: SingleCellExperiment 
# dim: 21369 109142 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(21369): ENSMMUG00000023296 ZNF692 ... ND6 CYTB
# rowData names(7): source type ... gene_biotype Symbol.uniq
# colnames(109142): 1_AAACGAAAGCGCCGTT-1 1_AAACGAACAGCCGTTG-1 ...
# 35_TTTGGTTGTGGCTAGA-1 35_TTTGTTGCACATGGTT-1
# colData names(14): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

sce.baboon <- logNormCounts(sce.baboon)
sce.baboon
# class: SingleCellExperiment 
# dim: 20842 57564 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(20842): PITHD1 ACTL8 ... ENSPANG00000037790 ENSPANG00000045064
# rowData names(7): source type ... gene_biotype Symbol.uniq
# colnames(57564): 1_AAACCCAAGACTGTTC-1 1_AAACCCAAGATCCAAA-1 ...
# 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# colData names(22): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):


# save new sce objects
save(sce.human,file=here(processed_dir,"sce_human_normalized.rda"))
save(sce.macaque,file=here(processed_dir,"sce_macaque_normalized.rda"))
save(sce.baboon,file=here(processed_dir,"sce_baboon_normalized.rda"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
# BiocSingular           1.16.0    2023-04-25 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# bluster                1.10.0    2023-04-25 [2] Bioconductor
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
# dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
# edgeR                  3.42.4    2023-05-31 [2] Bioconductor
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
# GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
# ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.3.1)
# ggplot2              * 3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.3.1)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# igraph                 1.5.1     2023-08-10 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# limma                  3.56.2    2023-06-04 [2] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# metapod                1.8.0     2023-04-25 [2] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# ScaledMatrix           1.8.1     2023-05-03 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater               * 1.28.0    2023-04-25 [2] Bioconductor
# scran                * 1.28.2    2023-07-23 [2] Bioconductor
# scuttle              * 1.10.2    2023-08-03 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                  0.6.3     
