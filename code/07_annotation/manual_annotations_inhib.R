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
library(dplyr)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(BiocParallel)
library(bluster)

## save directories
plot_dir = here("plots", "07_annotation", "seurat_subclustering_batch_correction")
processed_dir = here("processed-data","07_annotation")

# save sce
seurat.excit <- readRDS(here("processed-data", "07_annotation", "seurat_v4", "excit.integrated.rds"))
seurat.inhib <- readRDS(here("processed-data", "07_annotation", "seurat_v4", "inhib.integrated.rds"))
seurat.other <- readRDS(here("processed-data", "07_annotation", "seurat_v4", "other.integrated.rds"))



# ===== Plot standard seurat clusters next to species and subregion =======

# excitatory
pdf(here(plot_dir,"excit", "UMAP_new_integrations.pdf"), width=20, height=5)
p1 <- DimPlot(seurat.excit, reduction = "umap", group.by = c("ident")) +
    theme(legend.position = "none")

p2 <- DimPlot(seurat.excit, reduction = "umap", group.by = c("species"))
p3 <- DimPlot(seurat.excit, reduction = "umap", group.by = c("Subregion"))
print(p1+p2+p3)
dev.off()

# inhibitory
pdf(here(plot_dir,"inhib", "UMAP_new_integrations.pdf"), width=20, height=5)
p1 <- DimPlot(seurat.inhib, reduction = "umap", group.by = c("ident")) +
    theme(legend.position = "none")

p2 <- DimPlot(seurat.inhib, reduction = "umap", group.by = c("species"))
p3 <- DimPlot(seurat.inhib, reduction = "umap", group.by = c("Subregion"))
print(p1+p2+p3)
dev.off()

# other
pdf(here(plot_dir,"other", "UMAP_new_integrations.pdf"), width=20, height=5)
p1 <- DimPlot(seurat.other, reduction = "umap", group.by = c("ident")) +
    theme(legend.position = "none")

p2 <- DimPlot(seurat.other, reduction = "umap", group.by = c("species"))
p3 <- DimPlot(seurat.other, reduction = "umap", group.by = c("Subregion"))
print(p1+p2+p3)
dev.off()


# ========= Convert to SCE ==========

sce.inhib <- as.SingleCellExperiment(seurat.inhib, assay = "originalexp")
sce.inhib

sce.excit <- as.SingleCellExperiment(seurat.excit, assay = "originalexp")
sce.excit

sce.other <- as.SingleCellExperiment(seurat.other, assay = "originalexp")
sce.other

# ============ Inhibitroy Annotations ==============


# ===== marker genes =====
features <- c("GAD1", "GAD2",
              "SST", "PVALB",
              "LAMP5", "VIP",
              "NPY", "CORT",
              "CCK","CRH",
              "PDYN", "PENK",
              "PRKCD","CRHBP",
              "GAL", "DRD1", 
              "DRD2", "FOXP2",
              "TSHZ1", "NTS", 
              "NXPH1","CARTPT",
              "SNCG", "NOS1",
              "SLC17A7", "SLC17A6",
              "MBP", "GFAP",
              "PPP1R1B", "RELN",
              "TAC3", "CALB2",
              "ZFHX3"
              )

# drop features not in sce.inhib
features <- features[features %in% rownames(sce.inhib)]
features

pdf(here(plot_dir,"inhib", "Dotplot_inhib_marker_genes.pdf"), width=7.5, height=7.5)
plotDots(sce.inhib, 
         features = features, 
         group = "ident", 
         center=TRUE, 
         scale=TRUE)
dev.off()

# notes for cluster labels
# 0: LAMP5.1
# 1. TSHZ1.1
# 2. VIP.1
# 3. SST
# 4. VIP.2
# 5. PVALB
# 6. LAMP5.2
# 7. PENK_DRD2
# 8. Unknown.1
# 9. CCK
# 10. Unknown.2
# 11. TSHZ1.2
# 12. CARTPT.1
# 13. CARTPT.2
# 14. PVALB_2
# 15. PVALB.3
# 16. Unknown.3
# 17. PPP1R1B_DRD1
# 18. Unknown.5
# 19. SST_NOS1
# 20. Unknown.6

# make new $fine_type annotation column 
sce.inhib$fine_type <- c()

# assign annotations
sce.inhib$fine_type[sce.inhib$ident == 0] <- "LAMP5.1"
sce.inhib$fine_type[sce.inhib$ident == 1] <- "TSHZ1.1"
sce.inhib$fine_type[sce.inhib$ident == 2] <- "VIP.1"
sce.inhib$fine_type[sce.inhib$ident == 3] <- "SST"
sce.inhib$fine_type[sce.inhib$ident == 4] <- "VIP.2"
sce.inhib$fine_type[sce.inhib$ident == 5] <- "PVALB.1"
sce.inhib$fine_type[sce.inhib$ident == 6] <- "LAMP5.2"
sce.inhib$fine_type[sce.inhib$ident == 7] <- "PENK_DRD2"
sce.inhib$fine_type[sce.inhib$ident == 8] <- "Unknown.1"
sce.inhib$fine_type[sce.inhib$ident == 9] <- "CCK"
sce.inhib$fine_type[sce.inhib$ident == 10] <- "Unknown.2"
sce.inhib$fine_type[sce.inhib$ident == 11] <- "TSHZ1.2"
sce.inhib$fine_type[sce.inhib$ident == 12] <- "CARTPT.1"
sce.inhib$fine_type[sce.inhib$ident == 13] <- "CARTPT.2"
sce.inhib$fine_type[sce.inhib$ident == 14] <- "PVALB.2"
sce.inhib$fine_type[sce.inhib$ident == 15] <- "PVALB.3"
sce.inhib$fine_type[sce.inhib$ident == 16] <- "Unknown.3"
sce.inhib$fine_type[sce.inhib$ident == 17] <- "PPP1R1B_DRD1"
sce.inhib$fine_type[sce.inhib$ident == 18] <- "Unknown.5"
sce.inhib$fine_type[sce.inhib$ident == 19] <- "SST_NOS1"
sce.inhib$fine_type[sce.inhib$ident == 20] <- "Unknown.6"


sce.inhib$fine_type <- as.factor(sce.inhib$fine_type)

unique(sce.inhib$fine_type)

# remake umap with annotations
pdf(here(plot_dir,"inhib", "UMAP_inhib_annotated.pdf"), width = 10, height = 5)
plotReducedDim(sce.inhib, dimred="UMAP", colour_by = "fine_type", point_size=0.1, text_by="fine_type")
dev.off()


# ======== Dropping low quality cluster =====

# drop all cells in clusters 8, 13, 18, 20 from ident
dim(sce.inhib)
# [1] 13874 32698

sce.inhib <- sce.inhib[, sce.inhib$ident != 8]
sce.inhib <- sce.inhib[, sce.inhib$ident != 13]
sce.inhib <- sce.inhib[, sce.inhib$ident != 18]
sce.inhib <- sce.inhib[, sce.inhib$ident != 20]
dim(sce.inhib)
# [1] 13874 29361

# remake umap with annotations
pdf(here(plot_dir,"inhib", "UMAP_inhib_annotated_dropped.pdf"), width = 10, height = 5)
plotReducedDim(sce.inhib, dimred="UMAP", colour_by = "fine_type", point_size=0.1, text_by="fine_type")
dev.off()


# rerun UMAP

sce.inhib <- runUMAP(sce.inhib, dimred= "PCA", min_dist = 0.3)

pdf(here(plot_dir,"inhib", "NEW_UMAP_inhib_annotated_dropped.pdf"), width = 10, height = 5)
plotReducedDim(sce.inhib, dimred="UMAP", colour_by = "fine_type", point_size=0.1, text_by="fine_type")
dev.off()










