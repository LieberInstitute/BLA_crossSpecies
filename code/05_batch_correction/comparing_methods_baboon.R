# srun -n 1 --mem=100G --cpus-per-task=8 --pty bash -i
library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(batchelor)
library(harmony)
library(Seurat)
library(SeuratData) 

## save directories
plot_dir = here("plots", "05_batch_correction", "method_comparisons", 'Baboon')
processed_dir = here("processed-data","05_batch_correction")

# load sce
sce.baboon <- readRDS(file=here("processed-data","04_normalization", "sce.baboon_normalized.rds"))
sce.baboon
# class: SingleCellExperiment 
# dim: 13874 47039 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(5): ID Symbol Type gene_id_ncbi Symbol.uniq
# colnames(47039): 1_AAACCCAAGACTGTTC-1 1_AAACCCAAGATCCAAA-1 ...
# 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# colData names(22): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):


# baboon
set.seed(123) # for reproducibility
baboon.subset <- sce.baboon[, sample(ncol(sce.baboon), 20000)]
dim(baboon.subset)
# [1] 10973 20000

# ============== PCA and TSNE before correction ==================

# ===== baboon =====
set.seed(915)

# calculating nullResiduals without batch
baboon.subset <- scry::nullResiduals(baboon.subset, assay="counts", type="deviance")
baboon.subset <- scater::runPCA(baboon.subset, ncomponents = 50,
                      ntop = 5000,
                      exprs_values = "binomial_deviance_residuals",
                      scale = TRUE, name = "GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())

pdf(here(plot_dir, "PCA_baboon_GLM-PCA_uncorrected.pdf"))
plotReducedDim(baboon.subset, dimred="GLM-PCA", colour_by="Sample", ) +
    theme(legend.position = "none")
dev.off()


# calculating nullResiduals with batch
baboon.subset <- scry::nullResiduals(baboon.subset, assay="counts", type="deviance", batch='Sample')
baboon.subset <- scater::runPCA(baboon.subset, ncomponents = 50,
                               ntop = 5000,
                               exprs_values = "binomial_deviance_residuals",
                               scale = TRUE, name = "GLM-PCA_wBatch",
                               BSPARAM = BiocSingular::RandomParam())

pdf(here(plot_dir, "PCA_baboon_GLM-PCA_wBatch_uncorrected.pdf"))
plotReducedDim(baboon.subset, dimred="GLM-PCA_wBatch", colour_by="Sample", ) +
    theme(legend.position = "none")
dev.off()


# calculating UMAP
baboon.subset <- runUMAP(baboon.subset, dimred="GLM-PCA",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "UMAP_baboon_GLM-PCA_uncorrected.pdf"))
plotReducedDim(baboon.subset, dimred="UMAP", colour_by="Sample") +
    theme(legend.position = "none")
dev.off()

# calculate TSNE
baboon.subset <- runTSNE(baboon.subset, dimred="GLM-PCA",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "TSNE_baboon_GLM-PCA_uncorrected.pdf"))
plotReducedDim(baboon.subset, dimred="TSNE", colour_by="Sample") +
  theme(legend.position = "none")
dev.off()


# now normal PCA and UMAP
baboon.subset <- scater::runPCA(baboon.subset, ncomponents = 50,
                      ntop = 5000,
                      exprs_values = "logcounts",
                      scale = TRUE, name = "PCA",
                      BSPARAM = BiocSingular::RandomParam())

pdf(here(plot_dir, "PCA_baboon_PCA_uncorrected.pdf"))
plotReducedDim(baboon.subset, dimred="PCA", colour_by="Sample" )
dev.off()


# calculating UMAP
baboon.subset <- runUMAP(baboon.subset, dimred="PCA",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "UMAP_baboon_PCA_uncorrected.pdf"))
plotReducedDim(baboon.subset, dimred="UMAP", colour_by="Sample")
dev.off()

# calculate TSNE
baboon.subset <- runTSNE(baboon.subset, dimred="PCA",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "TSNE_baboon_PCA_uncorrected.pdf"))
plotReducedDim(baboon.subset, dimred="TSNE", colour_by="Sample")
dev.off()

# plot TSNE and color with doublet score
pdf(here(plot_dir, "TSNE_baboon_PCA_uncorrected_doublet_score.pdf"))
plotTSNE(baboon.subset, colour_by="doubletScore")
dev.off()
