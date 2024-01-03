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
plot_dir = here("plots", "05_batch_correction", "method_comparisons", 'Macaque')
processed_dir = here("processed-data","05_batch_correction")

sce.macaque <- readRDS(file=here("processed-data", "04_normalization", "sce.macaque_normalized.rds"))
sce.macaque
# class: SingleCellExperiment 
# dim: 13874 109142 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(7): source type ... gene_biotype Symbol.uniq
# colnames(109142): 1_AAACGAAAGCGCCGTT-1 1_AAACGAACAGCCGTTG-1 ...
# 35_TTTGGTTGTGGCTAGA-1 35_TTTGTTGCACATGGTT-1
# colData names(14): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):



# macaque
set.seed(123) # for reproducibility
macaque.subset <- sce.macaque[, sample(ncol(sce.macaque), 20000)]
dim(macaque.subset)
# [1] 10973 20000

# ============== PCA and TSNE before correction ==================

# ===== macaque =====
set.seed(915)

# calculating nullResiduals without batch
macaque.subset <- scry::nullResiduals(macaque.subset, assay="counts", type="deviance")
macaque.subset <- scater::runPCA(macaque.subset, ncomponents = 50,
                      ntop = 5000,
                      exprs_values = "binomial_deviance_residuals",
                      scale = TRUE, name = "GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())

pdf(here(plot_dir, "PCA_macaque_GLM-PCA_uncorrected.pdf"))
plotReducedDim(macaque.subset, dimred="GLM-PCA", colour_by="Sample", )
dev.off()


# calculating nullResiduals with batch
macaque.subset <- scry::nullResiduals(macaque.subset, assay="counts", type="deviance", batch='Sample')
macaque.subset <- scater::runPCA(macaque.subset, ncomponents = 50,
                               ntop = 5000,
                               exprs_values = "binomial_deviance_residuals",
                               scale = TRUE, name = "GLM-PCA_wBatch",
                               BSPARAM = BiocSingular::RandomParam())

pdf(here(plot_dir, "PCA_macaque_GLM-PCA_wBatch_uncorrected.pdf"))
plotReducedDim(macaque.subset, dimred="GLM-PCA_wBatch", colour_by="Sample", )
dev.off()


# calculating UMAP
macaque.subset <- runUMAP(macaque.subset, dimred="GLM-PCA",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "UMAP_macaque_GLM-PCA_uncorrected.pdf"))
plotReducedDim(macaque.subset, dimred="UMAP", colour_by="Sample")
dev.off()

# calculate TSNE
macaque.subset <- runTSNE(macaque.subset, dimred="GLM-PCA",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "TSNE_macaque_GLM-PCA_uncorrected.pdf"))
plotReducedDim(macaque.subset, dimred="TSNE", colour_by="Sample")
dev.off()


# now normal PCA and UMAP
macaque.subset <- scater::runPCA(macaque.subset, ncomponents = 50,
                      ntop = 5000,
                      exprs_values = "logcounts",
                      scale = TRUE, name = "PCA",
                      BSPARAM = BiocSingular::RandomParam())

pdf(here(plot_dir, "PCA_macaque_PCA_uncorrected.pdf"))
plotReducedDim(macaque.subset, dimred="PCA", colour_by="Sample", )
dev.off()


# calculating UMAP
macaque.subset <- runUMAP(macaque.subset, dimred="PCA",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "UMAP_macaque_PCA_uncorrected.pdf"))
plotReducedDim(macaque.subset, dimred="UMAP", colour_by="Sample")
dev.off()

# calculate TSNE
macaque.subset <- runTSNE(macaque.subset, dimred="PCA",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "TSNE_macaque_PCA_uncorrected.pdf"))
plotReducedDim(macaque.subset, dimred="TSNE", colour_by="Sample")
dev.off()

# plot TSNE and color with doublet score
pdf(here(plot_dir, "TSNE_macaque_PCA_uncorrected_doublet_score.pdf"))
plotTSNE(macaque.subset, colour_by="doubletScore")
dev.off()
