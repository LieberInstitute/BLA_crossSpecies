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
plot_dir = here("plots", "05_batch_correction", "method_comparisons", 'Human')
processed_dir = here("processed-data","05_batch_correction")

# load sce
sce.human <- readRDS(file=here("processed-data","04_normalization", "sce.human_normalized.rds"))
sce.human
# class: SingleCellExperiment 
# dim: 13874 21268 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(21268): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 5_TTTGTTGTCGGACTGC-1 5_TTTGTTGTCGTTGTTT-1
# colData names(28): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):


# ============== PCA and TSNE before correction ==================

# ===== HUMAN =====
set.seed(915)
human.subset <- sce.human

# calculating nullResiduals without batch
human.subset <- scry::nullResiduals(human.subset, assay="counts", type="deviance")
human.subset <- scater::runPCA(human.subset, ncomponents = 50,
                      ntop = 5000,
                      exprs_values = "binomial_deviance_residuals",
                      scale = TRUE, name = "GLM-PCA",
                      BSPARAM = BiocSingular::RandomParam())

pdf(here(plot_dir, "PCA_human_GLM-PCA_uncorrected.pdf"))
plotReducedDim(human.subset, dimred="GLM-PCA", colour_by="Sample", )
dev.off()


# calculating nullResiduals with batch
human.subset <- scry::nullResiduals(human.subset, assay="counts", type="deviance")
human.subset <- scater::runPCA(human.subset, ncomponents = 50,
                               ntop = 5000,
                               exprs_values = "binomial_deviance_residuals",
                               scale = TRUE, name = "GLM-PCA_wBatch",
                               BSPARAM = BiocSingular::RandomParam())

pdf(here(plot_dir, "PCA_human_GLM-PCA_wBatch_uncorrected.pdf"))
plotReducedDim(human.subset, dimred="GLM-PCA_wBatch", colour_by="Sample", )
dev.off()


# calculating UMAP
human.subset <- runUMAP(human.subset, dimred="GLM-PCA",
                        n_neighbors=5,
                        name="UMAP_5",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

human.subset <- runUMAP(human.subset, dimred="GLM-PCA",
                        n_neighbors=15,
                        name="UMAP_15",
                        BPPARAM = BiocParallel::MulticoreParam(workers=8))

human.subset <- runUMAP(human.subset, dimred="GLM-PCA",
                        n_neighbors=30,
                        name="UMAP_30",
                        BPPARAM = BiocParallel::MulticoreParam(workers=8))

human.subset <- runUMAP(human.subset, dimred="GLM-PCA",
                        n_neighbors=50,
                        name="UMAP_50",
                        BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "UMAP_human_GLM-PCA_neighbors_5.pdf"))
plotReducedDim(human.subset, dimred="UMAP_5", colour_by="Sample")
dev.off()

pdf(here(plot_dir, "UMAP_human_GLM-PCA_neighbors_15.pdf"))
plotReducedDim(human.subset, dimred="UMAP_15", colour_by="Sample")
dev.off()

pdf(here(plot_dir, "UMAP_human_GLM-PCA_neighbors_30.pdf"))
plotReducedDim(human.subset, dimred="UMAP_30", colour_by="Sample")
dev.off()

pdf(here(plot_dir, "UMAP_human_GLM-PCA_neighbors_50.pdf"))
plotReducedDim(human.subset, dimred="UMAP_50", colour_by="Sample")
dev.off()

# calculate TSNE
human.subset <- runTSNE(human.subset, dimred="GLM-PCA",
                        perplexity=2,
                        name="TSNE_2",
                        BPPARAM = BiocParallel::MulticoreParam(workers=8))

human.subset <- runTSNE(human.subset, dimred="GLM-PCA",
                        perplexity=5,
                        name="TSNE_5",
                        BPPARAM = BiocParallel::MulticoreParam(workers=8))

human.subset <- runTSNE(human.subset, dimred="GLM-PCA",
                        perplexity=10,
                        name="TSNE_10",
                        BPPARAM = BiocParallel::MulticoreParam(workers=8))

human.subset <- runTSNE(human.subset, dimred="GLM-PCA",
                        perplexity=25,
                        name="TSNE_25",
                        BPPARAM = BiocParallel::MulticoreParam(workers=8))

human.subset <- runTSNE(human.subset, dimred="GLM-PCA",
                        perplexity=40,
                        name="TSNE_40",
                        BPPARAM = BiocParallel::MulticoreParam(workers=8))

human.subset <- runTSNE(human.subset, dimred="GLM-PCA",
                        perplexity=500,
                        name="TSNE_500",
                        BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "TSNE_human_GLM-PCA_perplex_2.pdf"))
plotReducedDim(human.subset, dimred="TSNE_2", colour_by="Sample")
dev.off()

pdf(here(plot_dir, "TSNE_human_GLM-PCA_perplex_5.pdf"))
plotReducedDim(human.subset, dimred="TSNE_5", colour_by="Sample")
dev.off()

pdf(here(plot_dir, "TSNE_human_GLM-PCA_perplex_10.pdf"))
plotReducedDim(human.subset, dimred="TSNE_10", colour_by="Sample")
dev.off()

pdf(here(plot_dir, "TSNE_human_GLM-PCA_perplex_25.pdf"))
plotReducedDim(human.subset, dimred="TSNE_25", colour_by="Sample")
dev.off()

pdf(here(plot_dir, "TSNE_human_GLM-PCA_perplex_40.pdf"))
plotReducedDim(human.subset, dimred="TSNE_40", colour_by="Sample")
dev.off()

pdf(here(plot_dir, "TSNE_human_GLM-PCA_perplex_100.pdf"))
plotReducedDim(human.subset, dimred="TSNE_100", colour_by="Sample")
dev.off()

pdf(here(plot_dir, "TSNE_human_GLM-PCA_perplex_500.pdf"))
plotReducedDim(human.subset, dimred="TSNE_500", colour_by="Sample")
dev.off()




# now normal PCA and UMAP
human.subset <- scater::runPCA(human.subset, ncomponents = 50,
                      ntop = 5000,
                      exprs_values = "logcounts",
                      scale = TRUE, name = "PCA",
                      BSPARAM = BiocSingular::RandomParam())

pdf(here(plot_dir, "PCA_human_PCA_uncorrected.pdf"))
plotReducedDim(human.subset, dimred="PCA", colour_by="Sample", )
dev.off()


# calculating UMAP
human.subset <- runUMAP(human.subset, dimred="PCA",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "UMAP_human_PCA_uncorrected.pdf"))
plotReducedDim(human.subset, dimred="UMAP", colour_by="Sample")
dev.off()

# calculate TSNE
human.subset <- runTSNE(human.subset, dimred="PCA",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "TSNE_human_PCA_uncorrected.pdf"))
plotReducedDim(human.subset, dimred="TSNE", colour_by="Sample")
dev.off()

# plot TSNE and color with doublet score
pdf(here(plot_dir, "TSNE_human_PCA_uncorrected_doublet_score.pdf"))
plotTSNE(human.subset, colour_by="doubletScore")
dev.off()