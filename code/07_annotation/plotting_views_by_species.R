# srun -n 1 --mem=64G --cpus-per-task=1 --pty bash -i
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
library(cowplot)


## save directories
plot_dir = here("plots", "07_annotation", "final_views")
processed_dir = here("processed-data", "07_annotation")

# save sce
sce.excit <- readRDS(here(processed_dir, "sce.excit.integrated.annotated.rds"))
sce.inhib <- readRDS(here(processed_dir, "sce.inhib.final.rds"))
#sce.other <- readRDS(here(processed_dir, "sce.other.integrated.annotated.rds"))


# ========== Plotting views by species ===========

# === inhibitory ===
# subset to human, baboon, and macaque
inhib.human <- sce.inhib[,sce.inhib$species == "human"]
dim(inhib.human)
# [1] 13874  4355

inhib.baboon <- sce.inhib[,sce.inhib$species == "baboon"]
dim(inhib.baboon)
# [1] 13874  6301

inhib.macaque <- sce.inhib[,sce.inhib$species == "macaque"]
dim(inhib.macaque)
# [1] 13874 18705



# define consistent colors for subregions
mycolors <- c("Basal" = "#E41A1C", "Lateral" = "#377EB8", "Central Nucleus" = "#4DAF4A", "Accessory Basal" = "#984EA3", "NA" = "black")

# ploy
pdf(here(plot_dir, "inhibitory_views_clusters.pdf"), width=20, height=12)
p1 <- plotReducedDim(sce.inhib, dimred="UMAP", color_by="fine_type", point_size=0.1, text_by="fine_type") +
    theme_void() +
    ggtitle("All nuclei") +
    theme(plot.title = element_text(hjust = 0.5, size=22),
          legend.position="none")
p2 <- plotReducedDim(inhib.macaque, dimred="UMAP", color_by="fine_type", point_size=0.1,) +
    theme_void() +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=22)) + 
    ggtitle("Macaque")
p3 <- plotReducedDim(inhib.baboon, dimred="UMAP", color_by="fine_type", point_size=0.1,) +
    theme_void() +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=22)) +
    ggtitle("Baboon")
p4 <- plotReducedDim(inhib.human, dimred="UMAP", color_by="fine_type", point_size=0.1,) +
    theme_void() +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=22)) +
    ggtitle("Human")


p5 <- plotReducedDim(sce.inhib, dimred="UMAP", color_by="Subregion", point_size=0.1) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    theme(legend.position="none") 
p6 <- plotReducedDim(inhib.macaque, dimred="UMAP", color_by="Subregion", point_size=0.1) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    theme(legend.position="none") 
p7 <- plotReducedDim(inhib.baboon, dimred="UMAP", color_by="Subregion", point_size=0.1) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    theme(legend.position="none") 
p8 <- plotReducedDim(inhib.human, dimred="UMAP", color_by="Subregion", point_size=0.1) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    theme(legend.position="none") 

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol=4)
dev.off()



# ========= Excitatory cells ==========

# subset to human, baboon, and macaque
excit.human <- sce.excit[,sce.excit$species == "human"]
dim(excit.human)

excit.baboon <- sce.excit[,sce.excit$species == "baboon"]
dim(excit.baboon)

excit.macaque <- sce.excit[,sce.excit$species == "macaque"]
dim(excit.macaque)

# plot
pdf(here(plot_dir, "excitatory_views_by_species.pdf"), width=20, height=12)
p1 <- plotReducedDim(sce.excit, dimred="UMAP", color_by="fine_type", point_size=0.1, text_by="fine_type") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size=22),
          legend.position="none") +
    ggtitle("All nuclei")
p2 <- plotReducedDim(excit.macaque, dimred="UMAP", color_by="fine_type", point_size=0.1,) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size=22),
          legend.position="none") + 
    ggtitle("Macaque")
p3 <- plotReducedDim(excit.baboon, dimred="UMAP", color_by="fine_type", point_size=0.1,) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size=22),
          legend.position="none") +
    ggtitle("Baboon")
p4 <- plotReducedDim(excit.human, dimred="UMAP", color_by="fine_type", point_size=0.1,) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size=22),
          legend.position="none") +
    ggtitle("Human")

p5 <- plotReducedDim(sce.excit, dimred="UMAP", color_by="Subregion", point_size=0.1) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    theme(legend.position="none")
p6 <- plotReducedDim(excit.macaque, dimred="UMAP", color_by="Subregion", point_size=0.1) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    theme(legend.position="none")
p7 <- plotReducedDim(excit.baboon, dimred="UMAP", color_by="Subregion", point_size=0.1) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    theme(legend.position="none")
p8 <- plotReducedDim(excit.human, dimred="UMAP", color_by="Subregion", point_size=0.1) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    theme(legend.position="none")

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol=4)
dev.off()


