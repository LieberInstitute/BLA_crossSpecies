# srun -n 1 --mem=64G --cpus-per-task=8 --pty bash -i
remotes::install_version("Matrix", version = "1.6-1")
library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(batchelor)
library(harmony)
library(Seurat)


## save directories
plot_dir = here("plots", "05_batch_correction", "seurat_v4")
processed_dir = here("processed-data","05_batch_correction")

# load sce
# save combined, uncorrected sce
load(here("processed-data","04_norm_and_dim_reduction", "sce_combined_uncorrected.rda"))
combined

unique(combined$batch)
# [1] "Human"   "Macaque" "Baboon" 

combined$Species <- combined$batch
unique(combined$Species)


# ============ Batch correction with Seurat v4 (RPCA) ============
# remove duplicate columns 
combined<- combined[, !duplicated(colnames(combined))]
combined
# dim: 14391 187856 
# metadata(0):
#     assays(4): merged counts logcounts binomial_deviance_residuals
# rownames(14391): SAMD11 NOC2L ... EIF1AY RPS4Y2
# rowData names(4): gene_name gene_type Symbol.uniq gene_biotype
# colnames(187856): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# colData names(16): batch Sample ... sizeFactor Species
# reducedDimNames(4): PCA TSNE GLM-PCA_approx UMAP
# mainExpName: NULL
# altExpNames(0):

combined.seurat <- as.Seurat(combined, counts="counts", data="logcounts")
combined.seurat
# An object of class Seurat 
# 14391 features across 18000 samples within 1 assay 
# Active assay: originalexp (14391 features, 0 variable features)
# 2 layers present: counts, data
# 6 dimensional reductions calculated: PCA, TSNE, mnn_species_pca, mnn_samples_pca, mnn_species_sample_pca, mnn_samples_species_pca

# standard normalization in seurat
combined.seurat <- NormalizeData(combined.seurat)
combined.seurat <- FindVariableFeatures(combined.seurat)
combined.seurat <- ScaleData(combined.seurat)

# Run PCA
combined.seurat <- RunPCA(combined.seurat)

# Run UMAP
combined.seurat <- FindNeighbors(combined.seurat, dims=1:30, reduction="pca")
combined.seurat <- FindClusters(combined.seurat, resolution=2, cluster.name = "unintegrated_clusters")
combined.seurat <- RunUMAP(combined.seurat, dims=1:30, reduction="pca", reduction.name = "umap.unintegrated")

# plot uncorrected data
pdf(here(plot_dir, "UMAP_species_uncorrected.pdf"))
DimPlot(combined.seurat, reduction = "umap.unintegrated", group.by = c("Species"))
dev.off()


# split the dataset into a list of seurat objects by species
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

# find integration anchors
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features, reduction = "rpca")

# integrate the datasets
seurat.int <- IntegrateData(anchorset = seurat.anchors)
DefaultAssay(seurat.int) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.int <- ScaleData(seurat.int, verbose = FALSE)
seurat.int <- RunPCA(seurat.int, npcs = 30, verbose = FALSE)
seurat.int <- RunUMAP(seurat.int, reduction = "pca", dims = 1:30)
seurat.int <- FindNeighbors(seurat.int, reduction = "pca", dims = 1:30)
seurat.int <- FindClusters(seurat.int, resolution = 0.5)

# Visualization
pdf(here(plot_dir, "TSNE_species_integrated.pdf"))
DimPlot(seurat.int, reduction = "tsne", group.by = "Species")
dev.off()

pdf(here(plot_dir, "TSNE_clusters_integrated.pdf"))
DimPlot(seurat.int, reduction = "tsne",)
dev.off()


# save integrated seurat object
save(seurat.int, file = here(processed_dir, "seurat_integrated.rda"))
