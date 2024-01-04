# srun -n 1 --mem=400G --cpus-per-task=8 --pty bash -i
library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(batchelor)
library(harmony)
library(Seurat)
library(dplyr)


## save directories
plot_dir = here("plots", "05_batch_correction", "seurat_v4")
processed_dir = here("processed-data","05_batch_correction")

# load sce
# save combined, uncorrected sce
load(here("processed-data","05_batch_correction", "seurat_integrated.rda"))
seurat.int
# An object of class Seurat 
# 16391 features across 187856 samples within 2 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 2 layers present: data, scale.data
# 1 other assay present: originalexp
# 2 dimensional reductions calculated: pca, tsne

# change default assay to counts
DefaultAssay(seurat.int)
DefaultAssay(seurat.int) <- "originalexp"
DefaultAssay(seurat.int)


# ============ Run UMAP and plot initial Louvain clusters ============
# intial clustering was done during batch correction using Seurat's
# old Louvain algortihm (algorithm=1). Will additionally test the new Louvain 
# algorithm (algorithm=2) as well as the Leiden algorithm (algorithm=3). 

# split the dataset into a list of seurat objects by species
seurat.list <- SplitObject(seurat.int, split.by = "species")

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

# save integrated seurat object
save(seurat.anchors, file = here(processed_dir, "seurat.anchors_species.rda"))

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
pdf(here(plot_dir, "UMAP_species_integrated_by_species.pdf"))
DimPlot(seurat.int, reduction = "umap", group.by = "species")
dev.off()

pdf(here(plot_dir, "UMAP_samples_integrated_by_species.pdf"))
DimPlot(seurat.int, reduction = "umap", group.by = "Sample") + NoLegend()
dev.off()


# save integrated seurat object
save(seurat.int, file = here(processed_dir, "seurat_integrated_species.rda"))

# ======= high quality plotting =======

pdf(here(plot_dir, "UMAP_species_integrated_by_species.pdf"))
DimPlot(seurat.int, reduction = "umap", group.by = "species", raster=FALSE)
dev.off()

pdf(here(plot_dir, "UMAP_samples_integrated_by_species.pdf"))
DimPlot(seurat.int, reduction = "umap", group.by = "Sample", , raster=FALSE) + NoLegend()
dev.off()

