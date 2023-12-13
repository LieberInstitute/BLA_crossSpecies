# srun -n 1 --mem=64G --cpus-per-task=20 --pty bash -i
library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(batchelor)

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


which(seqnames(sce.macaque) == "MT")

rownames(sce.macaque)[which(seqnames(sce.macaque) == "MT")]

# use grep to get genes with "MT-" in their name in human
mt_human <- rownames(sce.human)[grep("MT-", rownames(sce.human))]
mt_human

# drop MT- from mt_human
mt_human <- gsub("MT-", "", mt_human)
mt_human


# use grep to get genes with "MT-" in their name in macaque
mt_macaque <- rownames(sce.macaque)[grep("MT-", rownames(sce.macaque))]
mt_macaque

# use grep to get genes with "MT-" in their name in baboon
mt_baboon <- rownames(sce.baboon)[grep("MT-", rownames(sce.baboon))]
mt_baboon

# =========== Calculate size factors ==============

human_v_mac <- intersect(rownames(sce.human), rownames(sce.macaque))
all_common_genes <- intersect(human_v_mac, rownames(sce.baboon))

sce.human <- sce.human[all_common_genes,]
sce.macaque <- sce.macaque[all_common_genes,]
sce.baboon <- sce.baboon[all_common_genes,]

dim(sce.human)
# [1] 20842 21268
dim(sce.macaque)
# [1]  14391 109142

# ===== Perform multi-experiment normalization =====

out <- multiBatchNorm(sces = list(sce.human, sce.macaque, sce.baboon))
sce.human <- out[[1]]
sce.macaque <- out[[2]]
sce.baboon <- out[[3]]


# ===== Feature seleection =====

dec1 <- modelGeneVar(sce.human)
dec2 <- modelGeneVar(sce.macaque)
dec3 <- modelGeneVar(sce.baboon)

combined.dec <- combineVar(dec1, dec2, dec3)
chosen.hvgs <- getTopHVGs(combined.dec, n=5000)

# As a sanity check, let's confirm that there are indeed batch effects in the data

combined <- correctExperiments(Human=sce.human, Macaque=sce.macaque, Baboon=sce.baboon, PARAM=NoCorrectParam())

set.seed(100)
combined <- runPCA(combined, subset_row=chosen.hvgs,
                   BPPARAM = BiocParallel::MulticoreParam(workers=20))
combined <- runTSNE(combined, dimred="PCA",
                    BPPARAM = BiocParallel::MulticoreParam(workers=20))

pdf(here(plot_dir, "TSNE_batch_effects.pdf"))
plotTSNE(combined, colour_by="batch")
dev.off()



# ============== Batch correction using MNN ============

set.seed(101)
# run with 20 cores
combined <- scry::nullResiduals(combined, assay="counts", type="deviance", )

combined <- runPCA(combined, ncomponents = 50,
       ntop = 3000,
       exprs_values = "binomial_deviance_residuals",
       scale = TRUE, name = "GLM-PCA_approx",
       BPPARAM = BiocParallel::MulticoreParam(workers=20))

pdf(here(plot_dir, "GLM_PCA_approx.pdf"))
plotReducedDim(combined, dimred = "GLM-PCA_approx", colour_by = "batch")
dev.off()



# calculate UMAP and TSNE using GLM-PCA_approx
set.seed(1234)
combined <- runUMAP(combined,
                    dimred = "GLM-PCA_approx",
                    name = "UMAP",
                    BPPARAM = BiocParallel::MulticoreParam(workers=20))

combined <- runTSNE(combined,
                    dimred = "GLM-PCA_approx",
                    name = "TSNE",
                    BPPARAM = BiocParallel::MulticoreParam(workers=20))



# plot UMAP and TSNE with color by species
pdf(here(plot_dir, "UMAP_uncorrected_species.pdf"))
plotUMAP(combined, colour_by = "batch", point_alpha = 0.2)
dev.off()

pdf(here(plot_dir, "TSNE_uncorrected_species.pdf"))
plotTSNE(combined, colour_by = "batch", point_alpha = 0.2)
dev.off()


# ploy UMAP and TSNE with color by library size
pdf(here(plot_dir, "UMAP_uncorrected_lib_size.pdf"))
plotUMAP(combined, colour_by = "sum", point_alpha = 0.2)
dev.off()

pdf(here(plot_dir, "TSNE_uncorrected_lib_size.pdf"))
plotTSNE(combined, colour_by = "sum", point_alpha = 0.2)
dev.off()


# ploy UMAP and TSNE with color by doublet score
pdf(here(plot_dir, "UMAP_uncorrected_doubletScore.pdf"))
plotUMAP(combined, colour_by = "doubletScore", point_alpha = 0.2)
dev.off()

pdf(here(plot_dir, "TSNE_uncorrected_doubletScore.pdf"))
plotTSNE(combined, colour_by = "doubletScore", point_alpha = 0.2)
dev.off()

combined

# save combined, uncorrected sce
save(combined,file = here(processed_dir,"sce_combined_uncorrected.rda"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()