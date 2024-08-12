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
plot_dir = here("plots", "07_annotation","06_excit_annotations")
#processed_dir = here("processed-data","07_annotation")
processed_dir <- here("processed-data")

# load sce
sce <- readRDS(here("processed-data", "sce_excit_final_subclusters_annotated.rds"))
sce

# ======= UMAPs ========


#subset to different species
sce_macaque <- sce[,which(colData(sce)$species == "macaque")]
sce_baboon <- sce[,which(colData(sce)$species == "baboon")]
sce_human <- sce[,which(colData(sce)$species == "human")]


dim(sce)
# [1] 13842 76519

dim(sce_macaque)
# [1] 13842 29206

dim(sce_baboon)
# [1] 13842  34833

dim(sce_human)
# [1] 13842 12480



# ====== UMAPs =======

# Define the number of colors you want
nb.cols <- length(unique(sce$fine_celltype))
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

png(here(plot_dir, "UMAP_excit_annotated_celltypes.png"), width=7, height=7, units="in", res=300)
p1 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "fine_celltype", text_by="fine_celltype", point_size=0.75) +
    scale_color_manual(values = mycolors) +
    theme_void() +
    theme(legend.position="none") 

p1
dev.off()



unique(sce_baboon$Subregion)
# [1] "Basal"   "Lateral"

unique(sce_macaque$Subregion)
# [1] "Lateral"         "Basal"           "Central Nucleus" "Accessory Basal"

# set colors so that they are consistent across species
subregion_colors <- c( "Accessory Basal" ="#8A2BE2", "Basal" = "#FF4500","Central Nucleus" = "#FFD700", "Lateral" = "#00BFFF")


# plot UMAP with idents
png(here(plot_dir, "UMAP_excit_across_species_subregions.png"), width=12, height=4, units="in", res=300)

p1 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = subregion_colors) +
    theme_void() +
    ggtitle("All nuclei (n = 76,519)")  +
    theme(legend.position="none")

p2 <- plotReducedDim(sce_macaque, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = subregion_colors) +
    theme_void() +
    ggtitle("Macaque (n = 29,206)")  +
    theme(legend.position="none")

p3 <- plotReducedDim(sce_baboon, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = subregion_colors) +
    theme_void() +
    ggtitle("Baboon (n = 34,833)") +
    theme(legend.position="none")

p4 <- plotReducedDim(sce_human, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = subregion_colors) +
    theme_void() +
    ggtitle("Human (n = 12,480)")  +
    theme(legend.position="none")


ggarrange(p1, p2, p3, p4, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()



# ========= Heatmap ==========

genes <- c(
    "ADARB2", "TRPS1", 
    "ESR1", "ADRA1A",
    "GRIK3", "TNS3", 
    "GULP1", "TRHDE", 
    "MEIS1", "PARD3B", 
    "MEIS2", "COL25A1", 
    "PEX5L", "MYRIP", 
    "RXFP1", "KIAA1217", 
    "SATB2", "MPPED1", 
    "SLC17A8", "ST8SIA2", 
    "ST18", "ABCA8", 
    "ZBTB20", "SLC4A4"
)


# == Create scaled matrix of average logocunts per celltype == 
out <- aggregateAcrossCells(sce, 
                            ids=sce$fine_celltype, 
                            subset_row=genes,
                            statistics="mean", 
                            use.assay.type="logcounts"
)

mat <- as.matrix(logcounts(out))
scaled_mat <- scale(t(mat))

# == Creating stacked bar plot species annotations == 
species_colors <- c("baboon" ="#fe9380", "human" = "#b0d5f5", "macaque" = "#ffc84d")



cell_types <- sce$fine_celltype
species <- sce$species

bar_width <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

species_proportions <- table(cell_types, species) / rowSums(table(cell_types, species))
ha_species <- rowAnnotation(Species = anno_barplot(species_proportions, gp = gpar(fill = species_colors),
                                                   bar_width=.8
), 
border=FALSE,
width=unit(2,"cm")
)



# == Creating stacked bar plot subregion annotations == 
sce.macaque <- sce[,which(colData(sce)$species == "macaque")]
cell_types <- sce.macaque$fine_celltype
subregions <- sce.macaque$Subregion  

bar_width <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

subregion_proportions <- table(cell_types, subregions) / rowSums(table(cell_types, subregions))
ha_subregion <- rowAnnotation(Subregion = anno_barplot(subregion_proportions, gp = gpar(fill = subregion_colors),
                                                       bar_width=.8
), 
border=FALSE,
width=unit(2,"cm"),
annotation_legend_param = list(title = "Subregion Proportions")
)


# == Creating celltype label annotations ===
# Create a named vector mapping cell types to colors

colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#FFA500", 
            "#800080", "#008000", "#000080", "#FFC0CB", "#D2B48C")

celltype_colors <- setNames(colors, unique(sce$fine_celltype))

anno_df <- data.frame(Cell_Type = out$fine_celltype)
ha_celltypes <- rowAnnotation(Cell_Type = anno_df$Cell_Type, 
                              col = list(Cell_Type = celltype_colors),
                              annotation_legend_param = list(Cell_Type = list(title = "Cell Types")),
                              show_legend=FALSE
)


# Custom legends
lgd_species <- Legend(labels = names(species_colors), title = "Species", 
                      legend_gp = gpar(fill = species_colors), pch = 16)

lgd_subregion <- Legend(labels = names(subregion_colors), title = "Subregion",
                        legend_gp = gpar(fill = subregion_colors), pch = 16)

lgd_list <- list(lgd_species, lgd_subregion)


#Heatmap Scale
color_scale <- colorRamp2(c(-4, 0, max(scaled_mat)), c("blue", "#EEEEEE", "red"))


ht <- Heatmap(scaled_mat,
              cluster_rows=FALSE,
              cluster_columns=FALSE,
              col=color_scale,
              rect_gp = gpar(col = "black", lwd = .75),
              left_annotation=ha_celltypes,
              show_heatmap_legend = FALSE,
              row_names_side="left"
) +
    ha_species +
    ha_subregion

ht_legend <- Legend(at = c(-4, 0, max(4)), col_fun = color_scale, title = "centered.scaled")


lgd_list <- packLegend(ht_legend, lgd_species, lgd_subregion, direction = "vertical")


png(here(plot_dir, "Heatmap_excitatory_celltypes.png"), width=10, height=6, units="in", res=300)
draw(ht, annotation_legend_list = lgd_list)
dev.off()










