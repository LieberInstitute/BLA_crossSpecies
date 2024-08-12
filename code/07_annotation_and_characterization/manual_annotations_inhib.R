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

sce.inhib$ident <- droplevels(sce.inhib$ident)


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
sce.inhib$fine_type[sce.inhib$ident == 9] <- "CCK"
sce.inhib$fine_type[sce.inhib$ident == 10] <- "ZFHX3"
sce.inhib$fine_type[sce.inhib$ident == 11] <- "TSHZ1.2"
sce.inhib$fine_type[sce.inhib$ident == 12] <- "CARTPT"
sce.inhib$fine_type[sce.inhib$ident == 14] <- "PVALB.2"
sce.inhib$fine_type[sce.inhib$ident == 15] <- "PVALB.3"
sce.inhib$fine_type[sce.inhib$ident == 16] <- "RELN"
sce.inhib$fine_type[sce.inhib$ident == 17] <- "PPP1R1B_DRD1"
sce.inhib$fine_type[sce.inhib$ident == 19] <- "SST_NOS1"


sce.inhib$fine_type <- as.factor(sce.inhib$fine_type)

unique(sce.inhib$fine_type)
# [1] LAMP5.2      VIP.2        LAMP5.1      TSHZ1.1      SST_NOS1     VIP.1        PVALB.1     
# [8] PENK_DRD2    SST          CCK          CARTPT       PVALB.2      PVALB.3      PPP1R1B_DRD1
# [15] ZFHX3        TSHZ1.2      RELN        
# 17 Levels: CARTPT CCK LAMP5.1 LAMP5.2 PENK_DRD2 PPP1R1B_DRD1 PVALB.1 PVALB.2 PVALB.3 RELN ... ZFHX3

# # remake umap with annotations
pdf(here(plot_dir,"inhib", "UMAP_inhib_annotated.pdf"), width = 10, height = 5)
plotReducedDim(sce.inhib, dimred="UMAP", colour_by = "fine_type", point_size=0.1, text_by="fine_type")
dev.off()



# rerun UMAP after dropping cells

set.seed(3500)
sce.inhib <- runUMAP(sce.inhib, dimred= "PCA", min_dist = 0.3)

pdf(here(plot_dir,"inhib", "NEW_UMAP_inhib_annotated_dropped.pdf"), width = 8, height = 5)
plotReducedDim(sce.inhib, dimred="UMAP", colour_by = "fine_type", point_size=0.1, text_by="fine_type")
dev.off()

pdf(here(plot_dir,"inhib", "UMAP_inhib_annotated_views.pdf"), width = 20, height = 5)
p1 <- plotReducedDim(sce.inhib, dimred="UMAP", colour_by = "fine_type", point_size=0.1, text_by="fine_type") +
  ggtitle("Clusters") +
  theme(legend.position = "none") 
p2 <- plotReducedDim(sce.inhib, dimred="UMAP", colour_by = "species", point_size=0.1) +
  ggtitle("Species")
p3 <- plotReducedDim(sce.inhib, dimred="UMAP", colour_by = "Subregion", point_size=0.1) +
  ggtitle("Subegion")

print(p1+p2+p3)
dev.off()

# save final sce inhib
saveRDS(sce.inhib, file = here(processed_dir, "sce.inhib.final.rds"))



# =========== Excitatory Annotations ==============


# ===== marker genes =====
features <- c("SLC17A7", "SLC17A6",
              "LAMP5","COL25A1",
              "GULP1","MOXD1",
              "RGS6", "CHGB",
              "RFX3"
)

# 0: GULP1_ALCAM
# 1: COL25A1_VWC2
# 2: GULP1_LAMP5
# 5: DAB1_MOXD1
# 6: SUSD4_GRIK1
# 8: GULP1_CPNE4
# 13: DGKG_TENM3
# 14: CUX2_TMEM132D
# 15: CHRM3_CPNE4

# features from above 

features <- c("LAMP5","COL25A1",
              "GULP1","MOXD1",
              "RGS6", "CHGB",
              "RFX3","ALCAM", 
              "VWC2", "GRIK1",
              "CPNE4", "TENM3",
              "TMEM132D", "CHRM3"
)

# drop features not in sce.inhib
features <- features[features %in% rownames(sce.excit)]
features

pdf(here(plot_dir,"excit", "Dotplot_excit_marker_genes_top.pdf"), width=7.5, height=7.5)
plotDots(sce.excit, 
         features = features, 
         group = "ident", 
         center=TRUE, 
         scale=TRUE)
dev.off()


# collapse these clusters together (e.g., 4 and 7 into 1:
# 1,4,7
# 2,3,10
# 21,0
# 6,9,17

sce.excit$ident[sce.excit$ident == 4] <- 1
sce.excit$ident[sce.excit$ident == 7] <- 1

sce.excit$ident[sce.excit$ident == 3] <- 2
sce.excit$ident[sce.excit$ident == 10] <- 2

sce.excit$ident[sce.excit$ident == 21] <- 0

sce.excit$ident[sce.excit$ident == 6] <- 6
sce.excit$ident[sce.excit$ident == 9] <- 6
sce.excit$ident[sce.excit$ident == 17] <- 6



# drop all cells in clusters 11, 12, 16, 18, 19, 20, 22, 23
dim(sce.excit)

sce.excit <- sce.excit[, sce.excit$ident != 11]
sce.excit <- sce.excit[, sce.excit$ident != 12]
sce.excit <- sce.excit[, sce.excit$ident != 16]
sce.excit <- sce.excit[, sce.excit$ident != 18]
sce.excit <- sce.excit[, sce.excit$ident != 19]
sce.excit <- sce.excit[, sce.excit$ident != 20]
sce.excit <- sce.excit[, sce.excit$ident != 22]
sce.excit <- sce.excit[, sce.excit$ident != 23]

dim(sce.excit)


# drop unused factor levels from ident
sce.excit$ident <- droplevels(sce.excit$ident)


set.seed(498)
sce.excit <- runUMAP(sce.excit, dimred= "PCA", min_dist = 0.3)

pdf(here(plot_dir,"excit", "UMAP_excit_annotated_dropped.pdf"), width = 8, height = 5)
plotReducedDim(sce.excit, dimred="UMAP", colour_by = "ident", point_size=0.1, text_by="ident")
dev.off()



# lets try finding marker genes using seurat...

seurat.excit <- as.Seurat(seurat.excit, counts = "counts", data = "logcounts")
Idents(seurat.excit) <- seurat.excit@meta.data$ident

seurat_markers <- Seurat::FindAllMarkers(seurat.excit, test.use = "wilcox", only.pos = TRUE, 
    min.pct = 0.25, logfc.threshold = 0.25
)

(top2 <- seurat_markers %>% dplyr::group_by(cluster) %>% 
        dplyr::top_n(n = 2, wt = avg_log2FC))

# drop any duplicated genes in top2
top2 <- top2[!duplicated(top2$gene),]

seurat.excit <- ScaleData(seurat.excit)

pdf(here(plot_dir,"excit", "Seurat_heatmap_marker_genes_top3.pdf"), width=7.5, height=7.5)
Seurat::DoHeatmap(seurat.excit, features = top2$gene) + NoLegend()
dev.off()

pdf(here(plot_dir,"excit", "Seurat_dotplot_marker_genes_top3.pdf"), width=7.5, height=7.5)
Seurat::DotPlot(seurat.excit, features = top2$gene) + RotatedAxis()
dev.off()

# new labels for sce.excit#ident clusters
# 0: GULP1_NTNG1
# 1: COL25A1_PEX5L
# 2: GULP1_ZBTB20
# 5: CHRM3_TRPS1
# 6: GPC5_ESR1
# 8: SAMD5_NTNG1
# 13: TENM3_SATB2
# 14: ZNF804B_TTN
# 15: GRM8_ADARB2

sce.temp <- sce.excit
sce.excit$fine_type <- c()

# rename clusters
sce.excit$fine_type[sce.excit$ident == 0] <- "GULP1_NTNG1"
sce.excit$fine_type[sce.excit$ident == 1] <- "COL25A1_PEX5L"
sce.excit$fine_type[sce.excit$ident == 2] <- "GULP1_ZBTB20"
sce.excit$fine_type[sce.excit$ident == 5] <- "CHRM3_TRPS1"
sce.excit$fine_type[sce.excit$ident == 6] <- "GPC5_ESR1"
sce.excit$fine_type[sce.excit$ident == 8] <- "SAMD5_NTNG1"
sce.excit$fine_type[sce.excit$ident == 13] <- "TENM3_SATB2"
sce.excit$fine_type[sce.excit$ident == 14] <- "ZNF804B_TTN"
sce.excit$fine_type[sce.excit$ident == 15] <- "GRM8_ADARB2"


pdf(here(plot_dir,"excit", "NEW_UMAP_excit_annotated_final.pdf"), width = 8, height = 5)
plotReducedDim(sce.excit, dimred="UMAP", colour_by = "fine_type", point_size=0.1, text_by="fine_type")
dev.off()

pdf(here(plot_dir,"excit", "UMAP_excit_annotated_final_views.pdf"), width = 20, height = 5)
p1 <- plotReducedDim(sce.excit, dimred="UMAP", colour_by = "fine_type", point_size=0.1, text_by="fine_type") +
    theme(legend.position = "none")
p2 <- plotReducedDim(sce.excit, dimred="UMAP", colour_by = "species", point_size=0.1, text_by="species")
p3 <- plotReducedDim(sce.excit, dimred="UMAP", colour_by = "Subregion", point_size=0.1, text_by="Subregion")

print(p1+p2+p3)
dev.off()


# save final sce.excit
saveRDS(sce.excit, here(processed_dir, "sce.excit.integrated.annotated.rds"))

# ========== Glia (other) annotation =============

pdf(here(plot_dir,"other", "UMAP_other_annotated_dropped.pdf"), width = 8, height = 5)
plotReducedDim(sce.other, dimred="UMAP", colour_by = "ident", point_size=0.1, text_by="ident")
dev.off()









    
    
    
    

