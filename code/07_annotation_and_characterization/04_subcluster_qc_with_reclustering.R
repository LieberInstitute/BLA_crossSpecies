
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
plot_dir = here("plots", "07_annotation","03_fine_annotations")
processed_dir = here("processed-data","07_annotation")


seurat.excit <- readRDS(here("processed-data", "07_annotation", "seurat_v4", "excit.integrated.rds"))
seurat.inhib <- readRDS(here("processed-data", "07_annotation", "seurat_v4", "inhib.integrated.rds"))


# ===== Plotting subclustered UMAPs =====
# excitatory
png(here(plot_dir,"UMAP_excit_views.png"), width=20, height=5, units="in", res=300)
p1 <- DimPlot(seurat.excit, reduction = "umap", group.by = c("ident")) +
    theme(legend.position = "none")

p2 <- DimPlot(seurat.excit, reduction = "umap", group.by = c("species"))
p3 <- DimPlot(seurat.excit, reduction = "umap", group.by = c("Subregion"))
print(p1+p2+p3)
dev.off()

# inhibitory
png(here(plot_dir, "UMAP_inhib_views.png"), width=20, height=5, units="in", res=300)
p1 <- DimPlot(seurat.inhib, reduction = "umap", group.by = c("ident")) +
    theme(legend.position = "none")

p2 <- DimPlot(seurat.inhib, reduction = "umap", group.by = c("species"))
p3 <- DimPlot(seurat.inhib, reduction = "umap", group.by = c("Subregion"))
print(p1+p2+p3)
dev.off()

# Plotting UMAPs with ident labels

png(here(plot_dir, "UMAP_excit_idents.png"), width=10, height=10, units="in", res=300)
p1 <- DimPlot(seurat.excit, reduction = "umap", group.by = c("ident")) 
LabelClusters(plot = p1, id = "ident")
dev.off()

png(here(plot_dir, "UMAP_inhib_idents.png"), width=10, height=10, units="in", res=300)
p1 <- DimPlot(seurat.inhib, reduction = "umap", group.by = c("ident"))
LabelClusters(plot = p1, id = "ident")
dev.off()


# drop low quality "splatter" clusters
excit_to_drop <- c("19") 
inhib_to_drop <- c("23","16","15","18","24","21","19")

# Subset the Seurat object to keep only the clusters you want
seurat.excit <- subset(seurat.excit, idents = setdiff(Idents(seurat.excit), excit_to_drop))
seurat.inhib <- subset(seurat.inhib, idents = setdiff(Idents(seurat.inhib), inhib_to_drop))

png(here(plot_dir, "UMAP_excit_idents_filtered.png"), width=10, height=10, units="in", res=300)
p1 <- DimPlot(seurat.excit, reduction = "umap", group.by = c("ident")) 
LabelClusters(plot = p1, id = "ident")
dev.off()

png(here(plot_dir, "UMAP_inhib_idents_filtered.png"), width=10, height=10, units="in", res=300)
p1 <- DimPlot(seurat.inhib, reduction = "umap", group.by = c("ident"))
LabelClusters(plot = p1, id = "ident")
dev.off()


# ===== rerun dim red and clustering ====


# excitatory
seurat.excit <- RunPCA(seurat.excit, npcs = 30, verbose = FALSE)
seurat.excit <- RunUMAP(seurat.excit, reduction = "pca", dims = 1:30)
seurat.excit <- FindNeighbors(seurat.excit, reduction = "pca", dims = 1:30)
seurat.excit <- FindClusters(seurat.excit, resolution = 0.4)

png(here(plot_dir, "UMAP_excit_idents_filtered_newUMAPs.png"), width=10, height=10, units="in", res=300)
p1 <- DimPlot(seurat.excit, reduction = "umap", group.by = c("ident")) 
LabelClusters(plot = p1, id = "ident")
dev.off()

# inhibitory
seurat.inhib <- RunPCA(seurat.inhib, npcs = 30, verbose = FALSE)
seurat.inhib <- RunUMAP(seurat.inhib, reduction = "pca", dims = 1:30)
seurat.inhib <- FindNeighbors(seurat.inhib, reduction = "pca", dims = 1:30)
seurat.inhib <- FindClusters(seurat.inhib, resolution = 0.4)

png(here(plot_dir, "UMAP_inhib_idents_filtered_newUMAPs.png"), width=10, height=10, units="in", res=300)
p1 <- DimPlot(seurat.inhib, reduction = "umap", group.by = c("ident"))
LabelClusters(plot = p1, id = "ident")
dev.off()



# ===== Convert to SCE =====

# == Excitatory ==
# Extract the counts and logcounts from the originalexp assay
counts <- GetAssayData(seurat.excit, assay = "originalexp", slot = "counts")
logcounts <- GetAssayData(seurat.excit, assay = "originalexp", slot = "data")

# Create a SingleCellExperiment object
sce.excit <- SingleCellExperiment(
  assays = list(counts = counts, logcounts = logcounts),
  colData = seurat.excit@meta.data
)

# Copy over the dim reductions
reducedDim(sce.excit, "PCA") <- Embeddings(seurat.excit, "pca")
reducedDim(sce.excit, "UMAP") <- Embeddings(seurat.excit, "umap")
sce.excit$subcluster_idents <- Idents(seurat.excit)

# == Inhibitory ==
# Extract the counts and logcounts from the originalexp assay
counts <- GetAssayData(seurat.inhib, assay = "originalexp", slot = "counts")
logcounts <- GetAssayData(seurat.inhib, assay = "originalexp", slot = "data")

# Create a SingleCellExperiment object
sce.inhib <- SingleCellExperiment(
  assays = list(counts = counts, logcounts = logcounts),
  colData = seurat.inhib@meta.data
)

# Copy over the dim reductions
reducedDim(sce.inhib, "PCA") <- Embeddings(seurat.inhib, "pca")
reducedDim(sce.inhib, "UMAP") <- Embeddings(seurat.inhib, "umap")
sce.inhib$subcluster_idents <- Idents(seurat.inhib)



# == Test plotting ===
# inhibitory
png(here(plot_dir, "UMAP_inhib_final_sce.png"), width=5, height=5, units="in", res=300)
p1 <- plotReducedDim(sce.inhib, dimred = "UMAP", colour_by = "subcluster_idents", text_by="subcluster_idents", point_size=0.2) +
    theme(legend.position = "none")
p1
dev.off()

# excitatory
png(here(plot_dir, "UMAP_excit_final_sce.png"), width=5, height=5, units="in", res=300)
p1 <- plotReducedDim(sce.excit, dimred = "UMAP", colour_by = "subcluster_idents", text_by="subcluster_idents", point_size=0.2) +
    theme(legend.position = "none")
p1
dev.off()

# save sce

saveRDS(sce.excit, file = here(processed_dir, "sce_excit_final_subclusters.rds"))
saveRDS(sce.inhib, file = here(processed_dir, "sce_inhib_final_subclusters.rds"))



png(here(plot_dir, "UMAP_excit_ncald.png"), width=10, height=5, units="in", res=300)
p1 <- plotReducedDim(sce.excit, dimred = "UMAP", colour_by = "DV_axis", text_by="subcluster_idents", point_size=0.2)

p2 <- plotReducedDim(sce.excit, dimred = "UMAP", colour_by = "Subregion", point_size=0.2)
p1+p2
dev.off()