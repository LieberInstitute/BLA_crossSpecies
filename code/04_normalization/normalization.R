library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(batchelor)
library(orthogene)

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

load(file=here("processed-data","03_quality_control","Baboon","PerCellQC","sce_post_qc_baboon.rda"))
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


# # ====== Subset to only 1:1 orthologs ======

# method <- "gprofiler"

# baboon_orthos <- orthogene::convert_orthologs(gene_df = counts(sce.baboon),
#                                         gene_input = "rownames", 
#                                         gene_output = "rownames", 
#                                         input_species = "macaque",
#                                         output_species = "human",
#                                         non121_strategy = "drop_both_species",
#                                         method = method) 
# # =========== REPORT SUMMARY ===========
# # Total genes dropped after convert_orthologs :
# #     5,487 / 21,091 (26%)
# # Total genes remaining after convert_orthologs :
# #     15,604 / 21,091 (74%)

# macaque_orthos <- orthogene::convert_orthologs(gene_df = counts(sce.macaque),
#                                               gene_input = "rownames", 
#                                               gene_output = "rownames", 
#                                               input_species = "baboon",
#                                               output_species = "human",
#                                               non121_strategy = "drop_both_species",
#                                               method = method) 
# # =========== REPORT SUMMARY ===========
# # Total genes dropped after convert_orthologs :
# #     6,221 / 21,369 (29%)
# # Total genes remaining after convert_orthologs :
# #     15,148 / 21,369 (71%)

# nhp_orthologs <- intersect(rownames(baboon_orthos), rownames(macaque_orthos))
# length(nhp_orthologs)
# # [1] 14025

# # get only orthologs that are in the human dataset
# valid.human.orthologs <- nhp_orthologs[nhp_orthologs %in% rownames(sce.human)]
# length(valid.human.orthologs)
# # [1] 13948

# # get only orthologs that are in the macaque dataset
# valid.macaque.orthologs <- nhp_orthologs[nhp_orthologs %in% rownames(sce.macaque)]
# length(valid.macaque.orthologs)
# # [1] 13952

# # get only orthologs that are in the baboon dataset
# valid.baboon.orthologs <- nhp_orthologs[nhp_orthologs %in% rownames(sce.baboon)]
# length(valid.baboon.orthologs)
# # [1] 13937

# # get only orthologs common to all datasets
# human_v_mac.orthologs <- intersect(valid.human.orthologs, valid.macaque.orthologs)
# final.121.orthologs <- intersect(human_v_mac.orthologs, valid.baboon.orthologs)
# length(final.121.orthologs)
# # [1] 13874

# # subset to only 1:1 orthologs
# sce.human.ortho <- sce.human[final.121.orthologs ,]
# sce.human.ortho
# # class: SingleCellExperiment 
# # dim: 13874 21268 
# # metadata(1): Samples
# # assays(1): counts
# # rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# # rowData names(7): source type ... gene_type Symbol.uniq
# # colnames(21268): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# # 5_TTTGTTGTCGGACTGC-1 5_TTTGTTGTCGTTGTTT-1
# # colData names(27): Sample Barcode ... discard_auto doubletScore
# # reducedDimNames(0):
# #     mainExpName: NULL
# # altExpNames(0):

# sce.macaque.ortho <- sce.macaque[final.121.orthologs,]
# sce.macaque.ortho
# # dim: 13874 109142 
# # metadata(1): Samples
# # assays(1): counts
# # rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# # rowData names(7): source type ... gene_biotype Symbol.uniq
# # colnames(109142): 1_AAACGAAAGCGCCGTT-1 1_AAACGAACAGCCGTTG-1 ...
# # 35_TTTGGTTGTGGCTAGA-1 35_TTTGTTGCACATGGTT-1
# # colData names(13): Sample Barcode ... discard_auto doubletScore
# # reducedDimNames(0):
# #     mainExpName: NULL
# # altExpNames(0):

# sce.baboon.ortho <- sce.baboon[final.121.orthologs,]
# sce.baboon.ortho
# # class: SingleCellExperiment 
# # dim: 13874 47039 
# # metadata(1): Samples
# # assays(1): counts
# # rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# # rowData names(5): ID Symbol Type gene_id_ncbi Symbol.uniq
# # colnames(47039): 1_AAACCCAAGACTGTTC-1 1_AAACCCAAGATCCAAA-1 ...
# # 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# # colData names(21): Sample Barcode ... discard_auto doubletScore
# # reducedDimNames(0):
# #     mainExpName: NULL
# # altExpNames(0):

# # rename variables
# sce.human <- sce.human.ortho
# sce.macaque <- sce.macaque.ortho
# sce.baboon <- sce.baboon.ortho


# ========= Drop lowly expressed genes ==========

# drops genes that don't have 1 UMI in at least 5% of cells
num_reads <- 1
num_cells <- 0.01*ncol(sce.human)
keep <- which(DelayedArray::rowSums(counts(sce.human) >= num_reads ) >= num_cells)
sce.human <- sce.human[keep,]
dim(sce.human)



num_reads <- 1
num_cells <- 0.01*ncol(sce.macaque)
keep <- which(DelayedArray::rowSums(counts(sce.macaque) >= num_reads ) >= num_cells)
sce.macaque <- sce.macaque[keep,]
dim(sce.macaque)


num_reads <- 1
num_cells <- 0.01*ncol(sce.baboon)
keep <- which(DelayedArray::rowSums(counts(sce.baboon) >= num_reads ) >= num_cells)
sce.baboon <- sce.baboon[keep,]
dim(sce.baboon)



# =========== Quick cluster, compute size factors, and normalize ==============

set.seed(100)

# normalize human data
clust.human <- quickCluster(sce.human)
sce.human <- computeSumFactors(sce.human, cluster=clust.human, min.mean=0.1)
sce.human <- logNormCounts(sce.human)

# normalize macaque data
clust.macaque <- quickCluster(sce.macaque)
sce.macaque <- computeSumFactors(sce.macaque, cluster=clust.macaque, min.mean=0.1)
sce.macaque <- logNormCounts(sce.macaque)

# normalize baboon data
clust.baboon <- quickCluster(sce.baboon)
sce.baboon <- computeSumFactors(sce.baboon, cluster=clust.baboon, min.mean=0.1)
sce.baboon <- logNormCounts(sce.baboon)

# save sce objects
saveRDS(sce.human, file = here(processed_dir, "sce.human_normalized.rds"))
saveRDS(sce.macaque, file = here(processed_dir, "sce.macaque_normalized.rds"))
saveRDS(sce.baboon, file = here(processed_dir, "sce.baboon_normalized.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

