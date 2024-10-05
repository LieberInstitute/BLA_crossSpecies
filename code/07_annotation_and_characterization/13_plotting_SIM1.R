library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(ggpubr)
library(RColorBrewer)
library(ComplexHeatmap)
library(scuttle)
library(circlize)

## save directories
plot_dir = here("plots", "07_annotation_and_characterization","04_excit_annotations")
processed_dir <- here("processed-data")

# load sce
sce.excit <- readRDS(here("processed-data", "07_annotation_and_characterization", "sce_excit_final_subclusters_annotated.rds"))
sce.inhib <- readRDS(here("processed-data","07_annotation_and_characterization", "sce_inhib_final_subclusters_annotated.rds"))
sce.other <- readRDS(here("processed-data","07_annotation_and_characterization", "sce_other_final_celltypes.rds"))


# load new combined sce
sce <- readRDS(here("processed-data","07_annotation_and_characterization","sce_FINAL_all_celltypes.rds"))
sce

sce.old <- readRDS(here("processed-data", "07_annotation_and_characterization", "sce_broad_annotations.rds"))
sce.old

colnames(colData(sce))
#  [1] "orig.ident"             "nCount_originalexp"     "nFeature_originalexp"  
#  [4] "Sample_num"             "Sample"                 "Species"               
#  [7] "Subject"                "Sex"                    "Region"                
# [10] "Subregion"              "DV_axis"                "PI.NeuN"               
# [13] "batch"                  "Barcode"                "sum"                   
# [16] "detected"               "subsets_Mito_sum"       "subsets_Mito_detected" 
# [19] "subsets_Mito_percent"   "total"                  "high_mito"             
# [22] "low_lib"                "low_genes"              "discard_auto"          
# [25] "doubletScore"           "sizeFactor"             "species"               
# [28] "unintegrated_clusters"  "seurat_clusters"        "integrated_snn_res.0.5"
# [31] "ident"                  "key"                    "broad_celltype"        
# [34] "fine_celltype"   


# ======= Add UMAPs to new SCE object ========

# 1. Find matching nuclei keys
matching_keys <- intersect(sce.old$key, sce$key)

# 2. Subset sce.old and sce to only include the matching nuclei
sce.old_matching <- sce.old[, sce.old$key %in% matching_keys]
sce_matching <- sce[, sce$key %in% matching_keys]

# 3. Reorder sce.old to match sce
sce.old_matching <- sce.old_matching[, match(sce$key, sce.old_matching$key, nomatch = 0)]

# 4. Transfer UMAP coordinates from sce.old to sce
umap_total <- reducedDims(sce.old_matching)$UMAP
reducedDims(sce)$UMAP_new <- umap_total


rm(sce.old_matching, sce_matching, sce.old)


# ======== Plotting whole UMAP =======
library(patchwork)



# inhibitory cell type colors
sce.inhib <- sce[,which(colData(sce)$broad_celltype == "Inhibitory")]
colors <- pals::cols25()[1:18]
inhibitory_colors <- setNames(colors, unique(sce.inhib$fine_celltype))

# excitatory cell type colors
sce.excit <- sce[,which(colData(sce)$broad_celltype == "Excitatory")]
nb.cols <- length(unique(sce.excit$fine_celltype))
colors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)
excitatory_colors <- setNames(colors, unique(sce.excit$fine_celltype))

# other cell type colors
sce.other <- sce[,which(colData(sce)$broad_celltype == "Non-neuronal")]
greys <- gray.colors(6)
other_colors <- setNames(greys, unique(sce.other$fine_celltype))

# combine all colors
celltype_colors <- c(inhibitory_colors, excitatory_colors, other_colors)

png(here(plot_dir, "UMAPs_SIM1.png"), width=15, height=7.5, units="in", res=300)

p1 <- plotReducedDim(sce, dimred = "UMAP_new", colour_by = "broad_celltype", point_size=1) +
    theme_void() +
    ggtitle("Broad cell types")  +
    theme(axis.line = element_blank(), # Remove the default axis lines
    axis.text = element_blank(), # Remove axis text
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(size = 18),
    legend.position = c(0.85, 0.9),
    legend.text = element_text(size = 12), # Increase legend text size
    legend.title = element_blank()
   ) +
   guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))


p2 <- plotReducedDim(sce, dimred = "UMAP_new", colour_by = "fine_celltype", point_size=1) +
    scale_color_manual(values = celltype_colors) +
    theme_void() +
    ggtitle("Fine cell types")  +
    theme(legend.position="none",
    axis.line = element_blank(), # Remove the default axis lines
    axis.text = element_blank(), # Remove axis text
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(size = 18)
    ) +
   guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) 

(p1+p2)
dev.off()



# ======= PLotting violins =========

sce_macaque <- sce.excit[,which(colData(sce.excit)$species == "macaque")]
sce_baboon <- sce.excit[,which(colData(sce.excit)$species == "baboon")]
sce_human <- sce.excit[,which(colData(sce.excit)$species == "human")]



png(here(plot_dir, "Violins_SIM1.png"), width=14, height=8, units="in", res=300)
p1 <- plotExpression(sce_macaque, x="fine_celltype", features=c('SLC17A6', "SIM1", "OTP"), colour_by="fine_celltype", ncol=1) +
    theme(legend.position="none") +
    ggtitle("Macaque") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- plotExpression(sce_baboon, x="fine_celltype", features=c('SLC17A6', "SIM1", "OTP"), colour_by="fine_celltype", ncol=1) +
    theme(legend.position="none") +
    ggtitle("Baboon") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- plotExpression(sce_human, x="fine_celltype", features=c('SLC17A6', "SIM1", "OTP"), colour_by="fine_celltype", ncol=1) +
    ggtitle("Human") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

(p1+p2+p3)
dev.off()
