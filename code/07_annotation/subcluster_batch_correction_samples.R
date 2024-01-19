# srun -n 1 --mem=200G --cpus-per-task=8 --pty bash -i
library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(batchelor)
library(harmony)
library(Seurat)
#library(future)


## save directories
plot_dir = here("plots", "07_annotation")
processed_dir = here("processed-data","07_annotation")


sce.excit <- readRDS(here("processed-data", "07_annotation", "sce_excit_subclustering.rds"))
sce.inhib <- readRDS(here("processed-data", "07_annotation", "sce_inhib_subclustering.rds"))
sce.other <- readRDS(here("processed-data", "07_annotation", "sce_other_subclustering.rds"))


# ========= Conver to Seurat v4 ==========

sce.objects <- list("sce.excit", "sce.inhib", "sce.other")
names <- c("excit", "inhib", "other")

for (i in 1:length(sce.objects)) {

    sce <- get(sce.objects[[3]])
    sce
    
    # drop duplicated colData 
    sce <- sce[,!duplicated(colnames(sce))]
    
    seurat <- as.Seurat(sce, counts="counts", data="logcounts")
    
    # standard normalization in seurat
    seurat <- NormalizeData(seurat)
    seurat <- FindVariableFeatures(seurat, nfeatures = 5000)
    seurat <- ScaleData(seurat)
    
    # Run PCA
    seurat <- RunPCA(seurat)
    
    # Run UMAP
    seurat <- FindNeighbors(seurat, dims=1:30, reduction="pca")
    seurat <- FindClusters(seurat, resolution=2, cluster.name = "unintegrated_clusters")
    seurat <- RunUMAP(seurat, dims=1:30, reduction="pca", reduction.name = "umap.unintegrated")
    
    # save seurat object
    saveRDS(seurat, file = paste0(processed_dir, "/seurat_v4/", names[i], "_unintegrated.rds"))
    
    # plot uncorrected data
    pdf(here(plot_dir, "seurat_subclustering_batch_correction", names[i], "UMAP_species_uncorrected.pdf"))
    DimPlot(seurat, reduction = "umap.unintegrated", group.by = c("species"))
    dev.off()
    
    pdf(here(plot_dir, "seurat_subclustering_batch_correction", names[i], "UMAP_samples_uncorrected.pdf"))
    DimPlot(seurat, reduction = "umap.unintegrated", group.by = c("Sample")) + NoLegend()
    dev.off()
    
    # split the dataset into a list of seurat objects by species
    seurat.list <- SplitObject(seurat, split.by = "species")
    
    # normalize and identify variable features for each dataset independently
    seurat.list<- lapply(X = seurat.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
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
    save(seurat.anchors, file = paste0(processed_dir, "/seurat_v4/", names[i], ".anchors.rda"))
    
    # integrate datasets
    seurat.int <- IntegrateData(anchorset = seurat.anchors, dims = 1:30)
    
    # standard seurat workflow for norm and vis
    seurat.int <- ScaleData(seurat.int, verbose = FALSE)
    seurat.int <- RunPCA(seurat.int, npcs = 30, verbose = FALSE)
    seurat.int <- RunUMAP(seurat.int, reduction = "pca", dims = 1:30)
    seurat.int <- FindNeighbors(seurat.int, reduction = "pca", dims = 1:30)
    seurat.int <- FindClusters(seurat.int, resolution = 0.5)
    
    # save integrated seurat object
    saveRDS(seurat.int, file = paste0(processed_dir, "/seurat_v4/", names[i], ".integrated.rds"))
    
    # plot integrated data
    pdf(here(plot_dir, "seurat_subclustering_batch_correction", names[i], "UMAP_species_integrated.pdf"))
    print(DimPlot(seurat.int, reduction = "umap", group.by = c("species")))
    dev.off()
    
    pdf(here(plot_dir, "seurat_subclustering_batch_correction", names[i], "UMAP_samples_integrated.pdf"))
    print(DimPlot(seurat.int, reduction = "umap", group.by = c("Sample")) + NoLegend())
    dev.off()

}



# ===== load for re-plotting =======

seurat.excit <- readRDS(here("processed-data", "07_annotation", "seurat_v4", "excit.integrated.rds"))
seurat.inhib <- readRDS(here("processed-data", "07_annotation", "seurat_v4", "inhib.integrated.rds"))
seurat.other <- readRDS(here("processed-data", "07_annotation", "seurat_v4", "other.integrated.rds"))


# plot integrated data
pdf(here(plot_dir, "seurat_subclustering_batch_correction", "excit", "UMAP_species_integrated.pdf"))
print(DimPlot(seurat.excit, reduction = "umap", group.by = c("species")))
dev.off()

pdf(here(plot_dir, "seurat_subclustering_batch_correction", "excit", "UMAP_samples_integrated.pdf"))
print(DimPlot(seurat.excit, reduction = "umap", group.by = c("Sample")) + NoLegend())
dev.off()


pdf(here(plot_dir, "seurat_subclustering_batch_correction", "inhib", "UMAP_species_integrated.pdf"))
print(DimPlot(seurat.inhib, reduction = "umap", group.by = c("species")))
dev.off()

pdf(here(plot_dir, "seurat_subclustering_batch_correction", "inhib", "UMAP_samples_integrated.pdf"))
print(DimPlot(seurat.inhib, reduction = "umap", group.by = c("Sample")) + NoLegend())
dev.off()


pdf(here(plot_dir, "seurat_subclustering_batch_correction", "other", "UMAP_species_integrated.pdf"))
print(DimPlot(seurat.other, reduction = "umap", group.by = c("species")))
dev.off()

pdf(here(plot_dir, "seurat_subclustering_batch_correction", "other", "UMAP_samples_integrated.pdf"))
print(DimPlot(seurat.other, reduction = "umap", group.by = c("Sample")) + NoLegend())
dev.off()





# ========================================================

i <- 3

sce <- get(sce.objects[[3]])
sce

# drop duplicated colData 
sce <- sce[,!duplicated(colnames(sce))]

seurat <- as.Seurat(sce, counts="counts", data="logcounts")

# standard normalization in seurat
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, nfeatures = 500)
seurat <- ScaleData(seurat)

# Run PCA
seurat <- RunPCA(seurat)

# Run UMAP
seurat <- FindNeighbors(seurat, dims=1:30, reduction="pca")
seurat <- FindClusters(seurat, resolution=2, cluster.name = "unintegrated_clusters")
seurat <- RunUMAP(seurat, dims=1:30, reduction="pca", reduction.name = "umap.unintegrated")

# save seurat object
#saveRDS(seurat, file = paste0(processed_dir, "/seurat_v4/", names[i], "_unintegrated.rds"))

# plot uncorrected data
# pdf(here(plot_dir, "seurat_subclustering_batch_correction", names[i], "UMAP_species_uncorrected.pdf"))
# DimPlot(seurat, reduction = "umap.unintegrated", group.by = c("species"))
# dev.off()
# 
# pdf(here(plot_dir, "seurat_subclustering_batch_correction", names[i], "UMAP_samples_uncorrected.pdf"))
# DimPlot(seurat, reduction = "umap.unintegrated", group.by = c("Sample")) + NoLegend()
# dev.off()

unique(seurat@meta.data$species)
# [1] "human"   "baboon"  "macaque"

# split the dataset into a list of seurat objects by species
seurat.list <- SplitObject(seurat, split.by = "species")

# normalize and identify variable features for each dataset independently
seurat.list<- lapply(X = seurat.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 500)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 500)
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# find integration anchors
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features, reduction = "rpca")

# save integrated seurat object
#save(seurat.anchors, file = paste0(processed_dir, "/seurat_v4/", names[i], ".anchors.rda"))

# integrate datasets
seurat.int <- IntegrateData(anchorset = seurat.anchors, dims = 1:30)

# standard seurat workflow for norm and vis
seurat.int <- ScaleData(seurat.int, verbose = FALSE)
seurat.int <- RunPCA(seurat.int, npcs = 30, verbose = FALSE)
seurat.int <- RunUMAP(seurat.int, reduction = "pca", dims = 1:30)
seurat.int <- FindNeighbors(seurat.int, reduction = "pca", dims = 1:30)
seurat.int <- FindClusters(seurat.int, resolution = 0.5)

# save integrated seurat object
#saveRDS(seurat.int, file = paste0(processed_dir, "/seurat_v4/", names[i], ".integrated.rds"))

# plot integrated data
pdf(here(plot_dir, "seurat_subclustering_batch_correction", names[i], "UMAP_species_integrated_500hvg.pdf"))
print(DimPlot(seurat.int, reduction = "umap", group.by = c("species")))
dev.off()

pdf(here(plot_dir, "seurat_subclustering_batch_correction", names[i], "UMAP_samples_integrated_500hvg.pdf"))
print(DimPlot(seurat.int, reduction = "umap", group.by = c("Sample")) + NoLegend())
dev.off()
