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
library(dplyr)x


## save directories
plot_dir = here("plots", "06_clustering", "Final")
processed_dir = here("processed-data","06_clustering")

# load sce
# save combined, uncorrected sce
load(here("processed-data","05_batch_correction", "seurat_integrated_species.rda"))
seurat.int
# An object of class Seurat 
# 16391 features across 187856 samples within 2 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 2 layers present: data, scale.data
# 1 other assay present: originalexp
# 2 dimensional reductions calculated: pca, tsne


# ============ Run UMAP and plot initial Louvain clusters ============
# intial clustering was done during batch correction using Seurat's
# old Louvain algortihm (algorithm=1). Will additionally test the new Louvain 
# algorithm (algorithm=2) as well as the Leiden algorithm (algorithm=3). 

seurat.int <- RunUMAP(seurat.int, reduction = "pca", dims = 1:30)

pdf(here(plot_dir, "UMAP_species_integrated_.pdf"))
DimPlot(seurat.int, reduction = "umap", group.by = "species")
dev.off()

pdf(here(plot_dir, "UMAP_clusters_integrated.pdf"))
DimPlot(seurat.int, reduction = "umap")
dev.off()

# ========= New Louvain algorithm =========

seurat.int <- FindNeighbors(seurat.int, reduction = "pca", dims = 1:30)
seurat.int <- FindClusters(seurat.int, resolution = 0.5, algorithm = 2)
seurat.int

pdf(here(plot_dir, "UMAP_clusters_integrated_new_louvain.pdf"))
DimPlot(seurat.int, reduction = "umap")
dev.off()

# ========= Leiden algorithm =========
seurat.int <- FindClusters(seurat.int, resolution = 0.5, algorithm = 3)
seurat.int

pdf(here(plot_dir, "UMAP_clusters_integrated_leiden.pdf"))
DimPlot(seurat.int, reduction = "umap")
dev.off()

# ========= Save Seurat object =========
save(seurat.int, file = here(processed_dir, "seurat_integrated_final.rda"))



seurat_obj <- seurat.int
# Calculate count of each species within each cluster
species_counts <- seurat_obj@meta.data %>%
    group_by(cluster = Idents(seurat_obj), species) %>%
    summarise(Count = n(), .groups = 'drop')

# Calculate total count for each cluster
total_counts <- species_counts %>%
    group_by(cluster) %>%
    summarise(Total = sum(Count), .groups = 'drop')

# Join to calculate the percentage
species_percentages <- species_counts %>%
    left_join(total_counts, by = "cluster") %>%
    mutate(Percentage = (Count / Total) * 100)

# Generate bar plots for each cluster
pdf(here(plot_dir, "species_percentages.pdf"))
ggplot(species_percentages, aes(x = species, y = Percentage, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~ cluster) +
    theme_minimal() +
    labs(x = "species", y = "Percentage", title = "Species Composition in Each Cluster")
dev.off()


# Calculate count of each species within each cluster
species_counts <- seurat_obj@meta.data %>%
    group_by(cluster = Idents(seurat_obj), species) %>%
    summarise(Count = n(), .groups = 'drop')

# Calculate total count for each cluster
total_counts <- species_counts %>%
    group_by(cluster) %>%
    summarise(Total = sum(Count), .groups = 'drop')

# Join to calculate the percentage
species_percentages <- species_counts %>%
    left_join(total_counts, by = "cluster") %>%
    mutate(Percentage = (Count / Total) * 100)

# Generate a stacked bar plot with percentages
pdf(here(plot_dir, "species_percentages_stacked.pdf"), width = 10, height = 3)
ggplot(species_percentages, aes(x = cluster, y = Percentage, fill = species)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_minimal() +
    labs(x = "Cluster", y = "Percentage", title = "Species Composition in Each Cluster (Stacked Bars)")
dev.off()


# plotting some basic genes
features <- c("SNAP25","GFAP","MOBP", "SLC17A7", "GAD2", "FOXP2", "SST", "VIP",'PVALB',"LAMP5")

# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
pdf(here(plot_dir, "DotPlot_genes.pdf"), width = 10, height = 7.5)
DotPlot(seurat.int, features = features) + RotatedAxis()
dev.off()

pdf(here(plot_dir, "FeaturePlot_genes.pdf"), width = 20, height = 20)
FeaturePlot(seurat.int, features = features)
dev.off()