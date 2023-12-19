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
#library(SeuratWrappers)

## save directories
plot_dir = here("plots", "05_batch_correction", "method_comparisons")
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


# ===== Feature seleection =====

# # drops genes that don't have 1 UMI in at least 5% of cells
# num_reads <- 1
# num_cells <- 0.01*ncol(sce.human)
# keep <- which(DelayedArray::rowSums(counts(sce.human) >= num_reads ) >= num_cells)
# sce.human <- sce.human[keep,]
# dim(sce.human)
# # [1] 11779 21268 
# 
# 
# num_reads <- 1
# num_cells <- 0.01*ncol(sce.macaque)
# keep <- which(DelayedArray::rowSums(counts(sce.macaque) >= num_reads ) >= num_cells)
# sce.macaque <- sce.macaque[keep,]
# dim(sce.macaque)
# # [1] 10973 47039
# 
# num_reads <- 1
# num_cells <- 0.01*ncol(sce.baboon)
# keep <- which(DelayedArray::rowSums(counts(sce.baboon) >= num_reads ) >= num_cells)
# sce.baboon <- sce.baboon[keep,]
# dim(sce.baboon)
# # [1] 10973 47039



# ===== Subsetting to 10k cells =====

# # Initialize an empty list to store subsets
# subsets <- list()
# 
# # Loop through each batch
# for (batch_id in unique(combined$batch)) {
#     # Filter cells from the current batch
#     batch_cells <- combined[, combined$batch == batch_id]
#     
#     # Randomly sample 6000 nuclei from this batch
#     set.seed(123) # for reproducibility
#     sampled_cells <- batch_cells[, sample(ncol(batch_cells), 10000)]
#     
#     # Add the sampled cells to the list
#     subsets[[batch_id]] <- sampled_cells
# }

randomly subset human, baboon, and macaque data to 10000 cells

# human
set.seed(123) # for reproducibility
human.subset <- sce.human[, sample(ncol(sce.human), 20000)]
dim(human.subset)
# [1] 11779 10000

# macaque
set.seed(123) # for reproducibility
macaque.subset <- sce.macaque[, sample(ncol(sce.macaque), 20000)]
dim(macaque.subset)
# [1] 11173 10000

# baboon
set.seed(123) # for reproducibility
baboon.subset <- sce.baboon[, sample(ncol(sce.baboon), 20000)]
dim(baboon.subset)
# [1] 10973 10000



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
human.subset <- scry::nullResiduals(human.subset, assay="counts", type="deviance", batch='Sample')
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
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "UMAP_human_GLM-PCA_uncorrected.pdf"))
plotReducedDim(human.subset, dimred="UMAP", colour_by="Sample")
dev.off()

# calculate TSNE
human.subset <- runTSNE(human.subset, dimred="GLM-PCA",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "TSNE_human_GLM-PCA_uncorrected.pdf"))
plotReducedDim(human.subset, dimred="TSNE", colour_by="Sample")
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




# =========== Baboons =============

# calculating nullResiduals without batch
baboon.subset <- scry::nullResiduals(baboon.subset, assay="counts", type="deviance")
baboon.subset <- scater::runPCA(baboon.subset, ncomponents = 50,
                               ntop = 5000,
                               exprs_values = "binomial_deviance_residuals",
                               scale = TRUE, name = "GLM-PCA",
                               BSPARAM = BiocSingular::RandomParam())

pdf(here(plot_dir, "PCA_baboon_GLM-PCA_uncorrected.pdf"))
plotReducedDim(baboon.subset, dimred="GLM-PCA", colour_by="Sample", )
dev.off()

# ============== MNN on Species > Samples ============




# Run batch correction with MNN across Samples
mnn_species_pca <- batchelor::reducedMNN(reducedDim(combined.subset, "PCA"),
                                         batch=as.factor(combined.subset$Species))

reducedDim(combined.subset,"mnn_species_pca") <- mnn_species_pca$corrected


# Run batch correction with MNN across Species 
mnn_species_sample_pca <- batchelor::reducedMNN(reducedDim(combined.subset, "mnn_species_pca"),
                                                batch=as.factor(combined.subset$Sample))

reducedDim(combined.subset,"mnn_species_sample_pca") <- mnn_species_sample_pca$corrected


# compute TSNE
combined.subset <- runTSNE(combined.subset, dimred="mnn_species_sample_pca",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))


# plot
pdf(here(plot_dir, "TSNE_mnn_species_sample.pdf"))
plotReducedDim(combined.subset, dim="TSNE", colour_by="Species")
dev.off()



# ================ Hrmony on Species > Samples ==============

combined.subset <- RunHarmony(combined.subset, c("Species", "Sample"))
combined.subset

# compute TSNE
combined.subset <- runTSNE(combined.subset, dimred="HARMONY",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))


# plot
pdf(here(plot_dir, "TSNE_harmony_species_samples.pdf"))
plotReducedDim(combined.subset, dim="TSNE", colour_by="Species")
dev.off()


# ================ Hrmony on Samples > Species ==============

combined.subset <- RunHarmony(combined.subset, c("Species"))
combined.subset

# compute TSNE
combined.subset <- runTSNE(combined.subset, dimred="HARMONY",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))


# plot
pdf(here(plot_dir, "TSNE_harmony_species.pdf"))
plotReducedDim(combined.subset, dim="TSNE", colour_by="Species")
dev.off()


# ================ Seurat ==============
# remove duplicate columns 
combined.subset <- combined.subset[, !duplicated(colnames(combined.subset))]

combined.seurat <- as.Seurat(combined.subset, counts="counts", data="logcounts")
combined.seurat
# An object of class Seurat 
# 14391 features across 18000 samples within 1 assay 
# Active assay: originalexp (14391 features, 0 variable features)
# 2 layers present: counts, data
# 6 dimensional reductions calculated: PCA, TSNE, mnn_species_pca, mnn_samples_pca, mnn_species_sample_pca, mnn_samples_species_pca

DimPlot(combined.seurat, reduction="TSNE", group.by="Species")

combined.seurat <- NormalizeData(combined.seurat)
combined.seurat <- FindVariableFeatures(combined.seurat)
combined.seurat <- ScaleData(combined.seurat)

#combined.seurst <- SCTransform(combined.seurat, verbose = FALSE)
combined.seurat <- RunPCA(combined.seurat)

combined.seurat <- FindNeighbors(combined.seurat, dims=1:30, reduction="pca")
combined.seurat <- FindClusters(combined.seurat, resolution=2, cluster.name = "unintegrated_clusters")
combined.seurat <- RunUMAP(combined.seurat, dims=1:30, reduction="pca", reduction.name = "umap.unintegrated")


DimPlot(combined.seurat, reduction = "umap.unintegrated", group.by = c("Species"))


# split the dataset into a list of two seurat objects (stim and CTRL)
seurat.list <- SplitObject(combined.seurat, split.by = "Species")

# normalize and identify variable features for each dataset independently
seurat.list<- lapply(X = seurat.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = seurat.list)
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features, reduction = "rpca")

# this command creates an 'integrated' data assay
seurat.int <- IntegrateData(anchorset = seurat.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seurat.int) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.int <- ScaleData(seurat.int, verbose = FALSE)
seurat.int <- RunPCA(seurat.int, npcs = 30, verbose = FALSE)
seurat.int <- RunTSNE(seurat.int, reduction = "pca", dims = 1:30)
seurat.int <- FindNeighbors(seurat.int, reduction = "pca", dims = 1:30)
seurat.int <- FindClusters(seurat.int, resolution = 0.5)

# Visualization
pdf(here(plot_dir, "Seurat_UMAP_unintegrated.pdf"))
DimPlot(seurat.int, reduction = "tsne", group.by = "Species")
dev.off()


