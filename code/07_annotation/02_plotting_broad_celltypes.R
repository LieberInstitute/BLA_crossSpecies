# srun -n 1 --mem=32G --cpus-per-task=1 --pty bash -i
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
library(forcats)

## save directories
plot_dir = here("plots", "07_annotation", "01_broad_annotations")

# load sce
sce <- readRDS(here("processed-data", "07_annotation", "sce_broad_annotations.rds"))
sce

# ========== Comparing macaque, baboon, and human cells ========

#subset to different species
# sce_macaque <- sce[,which(colData(sce)$species == "macaque")]
# sce_baboon <- sce[,which(colData(sce)$species == "baboon")]
# sce_human <- sce[,which(colData(sce)$species == "human")]

# dim(sce)
# # [1]  13842 171928

# dim(sce_macaque)
# # [1]  13842 104797

# dim(sce_baboon)
# # [1] 13842 46111

# dim(sce_human)
# # [1] 13842 21020

png(here(plot_dir, "UMAP_across_species.png"), width=12, height=4, units="in", res=300)

p1 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "broad_celltype", point_size=0.1) +
    theme_void() +
    ggtitle("All nuclei (n = 171,928)")  +
    theme(legend.position="none") 

p2 <- plotReducedDim(sce_macaque, dimred = "UMAP", colour_by = "broad_celltype", point_size=0.1) +
    theme_void() +
    ggtitle("Macaque (n = 104,797)")  +
    theme(legend.position="none") 

p3 <- plotReducedDim(sce_baboon, dimred = "UMAP", colour_by = "broad_celltype", point_size=0.1) +
    theme_void() +
    ggtitle("Baboon (n = 46,111)") +
    theme(legend.position="none") 

p4 <- plotReducedDim(sce_human, dimred = "UMAP", colour_by = "broad_celltype", point_size=0.1) +
    theme_void() +
    ggtitle("Human (n = 21,020)")  +
    theme(legend.position="none") 


ggarrange(p1, p2, p3, p4, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()



# ======== Comparing macaque, baboon, and human cells across $ident =========


# Define the number of colors you want
nb.cols <- length(unique(sce$ident))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

# plot UMAP with idents
png(here(plot_dir, "UMAP_across_species_idents.png"), width=12, height=4, units="in", res=300)

p1 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "ident", point_size=0.1) +
    scale_fill_manual(values = mycolors) +
    theme_void() +
    ggtitle("All nuclei (n = 171,928)")  +
    theme(legend.position="none")

p2 <- plotReducedDim(sce_macaque, dimred = "UMAP", colour_by = "ident", point_size=0.1) +
    scale_fill_manual(values = mycolors) +
    theme_void() +
    ggtitle("Macaque (n = 104,797)")  +
    theme(legend.position="none")

p3 <- plotReducedDim(sce_baboon, dimred = "UMAP", colour_by = "ident", point_size=0.1) +
    scale_fill_manual(values = mycolors) +
    theme_void() +
    ggtitle("Baboon (n = 46,111)") +
    theme(legend.position="none")

p4 <- plotReducedDim(sce_human, dimred = "UMAP", colour_by = "ident", point_size=0.1) +
    scale_fill_manual(values = mycolors) +
    theme_void() +
    ggtitle("Human (n = 21,020)")  +
    theme(legend.position="none")


ggarrange(p1, p2, p3, p4, ncol=4, nrow=1)
dev.off()




# ========= Plotting UMAPs with subregions =========

unique(sce_baboon$Subregion)
# [1] "Basal"   "Lateral"

unique(sce_macaque$Subregion)
# [1] "Lateral"         "Basal"           "Central Nucleus" "Accessory Basal"

# set colors so that they are consistent across species
mycolors <- c("Basal" = "#E41A1C", "Lateral" = "#377EB8", "Central Nucleus" = "#4DAF4A", "Accessory Basal" = "#984EA3", "NA" = "black")

# plot UMAP with idents
png(here(plot_dir, "UMAP_across_species_subregions.png"), width=12, height=4, units="in", res=300)

p1 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    ggtitle("All nuclei (n = 171,928)")  +
    theme(legend.position="none")

p2 <- plotReducedDim(sce_macaque, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    ggtitle("Macaque (n = 104,797)")  +
    theme(legend.position="none")

p3 <- plotReducedDim(sce_baboon, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    ggtitle("Baboon (n = 46,111)") +
    theme(legend.position="none")

p4 <- plotReducedDim(sce_human, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    ggtitle("Human (n = 21,020)")  +
    theme(legend.position="none")


ggarrange(p1, p2, p3, p4, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()



# ========== Plot for Figure 1 ============

nb.cols <- length(unique(sce$ident))
ident_colors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

region_colors <- c("Basal" = "#E41A1C", "Lateral" = "#377EB8", "Central Nucleus" = "#4DAF4A", "Accessory Basal" = "#984EA3", "NA" = "black")

species_colors <- c("human" = "#375E97", "baboon" = "#FB6542", "macaque" = "#6AB187")

png(here(plot_dir, "UMAP_Figure1.png"), width=10, height=10, units="in", res=300)

# reorder to best show mixing across species
sce$species <- fct_relevel(sce$species, "macaque", "baboon", "human")
sce <- sce[,order(colData(sce)$species)]


p1 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "species", point_size=0.1) +
    scale_color_manual(values = species_colors) +
    theme_void() +
    ggtitle("Species")  +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_blank(), # Remove the default axis lines
    axis.text = element_blank(), # Remove axis text
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(size = 24),
    legend.position = c(0.85, 0.9),
    legend.text = element_text(size = 12), # Increase legend text size
    legend.title = element_blank()
   ) +
   guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))


# reorder again to best show differences across subregions
sce$species <- fct_relevel(sce$species, "human", "baboon", "macaque")
sce <- sce[,order(colData(sce)$species)]

p2 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = region_colors) +
    theme_void() +
    ggtitle("Amygdala subregion")  +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_blank(), # Remove the default axis lines
    axis.text = element_blank(), # Remove axis text
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(size = 24),
    legend.position = c(0.82, 0.85),
    legend.text = element_text(size = 12), # Increase legend text size
    legend.title = element_blank()
   ) +
   guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))


p3 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "broad_celltype", point_size=0.1) +
    theme_void() +
    ggtitle("Broad cell type clusters")  +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_blank(), # Remove the default axis lines
    axis.text = element_blank(), # Remove axis text
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(size = 24),
    legend.position = c(0.85, 0.9),
    legend.text = element_text(size = 12), # Increase legend text size
    legend.title = element_blank()
   ) +
   guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))


p4 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "ident", point_size=0.1) +
    scale_color_manual(values = ident_colors) +
    theme_void() +
    ggtitle("Fine cell type clusters")  +
    theme(legend.position="none",
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line = element_blank(), # Remove the default axis lines
    axis.text = element_blank(), # Remove axis text
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(size = 24)
    ) 

(p1+p2)/(p3+p4)
dev.off()
