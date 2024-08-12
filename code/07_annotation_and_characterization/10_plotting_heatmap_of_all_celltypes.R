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
plot_dir = here("plots", "07_annotation","broad_annotations")
#processed_dir = here("processed-data","07_annotation")
processed_dir <- here("processed-data")

# load sce
sce.excit <- readRDS(here("processed-data", "sce_excit_final_subclusters_annotated.rds"))
sce.inhib <- readRDS(here("processed-data", "sce_inhib_final_subclusters_annotated.rds"))
sce.other <- readRDS(here("processed-data", "JHPCE","sce_other_final_celltypes.rds"))

# combine SCE objects
sce.inhib$subcluster_idents <- NULL
sce.inhib$integrated_snn_res.0.4 <- NULL
sce.excit$integrated_snn_res.0.4 <- NULL

sce.neurons <- cbind(sce.excit, sce.inhib)
sce <- cbind(sce.neurons, sce.other)

rm(sce.excit, sce.inhib, sce.other, sce.neurons)
gc()



# ===== Dendrogram based on HVG expression =====

# Calculate HVGs
dec <- modelGeneVar(sce)
HVGs <- getTopHVGs(dec, n=500)
hvg_matrix <- logcounts(sce)[HVGs, ]

dist_matrix <- proxyC::dist(t(hvg_matrix), method = "euclidean") # Transpose because we want distance between cells
hvg_dend <- as.dendrogram(hclust(dist_matrix))



# ========= Heatmap ==========

genes <- c(
    #"SNAP25", 
    "SYT1",
    "SLC17A7",
    #"GAD1"
    "GAD2"
    #"SLC1A2", # astrocyte
    #"FOXJ1", # ependymal
    #"MBP", # oligodendrocyte
    #"TMEM119",  # microglia
    #"RGS5", # Endothelial
    #"CD9" # OPC
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
subregion_colors <- c( "Accessory Basal" ="#8A2BE2", "Basal" = "#FF4500","Central Nucleus" = "#FFD700", "Lateral" = "#00BFFF")




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
# 
# colors <- colorRampPalette(brewer.pal(36, "Set1"))(nb.cols)
# 
# celltype_colors <- setNames(colors, unique(sce$fine_celltype))
# 
# anno_df <- data.frame(Cell_Type = out$fine_celltype)
# ha_celltypes <- rowAnnotation(Cell_Type = anno_df$Cell_Type, 
#                               col = list(Cell_Type = celltype_colors),
#                               annotation_legend_param = list(Cell_Type = list(title = "Cell Types")),
#                               show_legend=FALSE
# )


# Custom legends
lgd_species <- Legend(labels = names(species_colors), title = "Species", 
                      legend_gp = gpar(fill = species_colors), pch = 16)

lgd_subregion <- Legend(labels = names(subregion_colors), title = "Subregion",
                        legend_gp = gpar(fill = subregion_colors), pch = 16)

lgd_list <- list(lgd_species, lgd_subregion)


#Heatmap Scale
color_scale <- colorRamp2(c(-4, 0, max(scaled_mat)), c("blue", "#EEEEEE", "red"))


ht <- Heatmap(scaled_mat,
              cluster_rows=TRUE,
              cluster_columns=FALSE,
              col=color_scale,
              rect_gp = gpar(col = "black", lwd = .75),
              #left_annotation=ha_celltypes,
              show_heatmap_legend = FALSE,
              row_names_side="left"
) +
    ha_species +
    ha_subregion

ht_legend <- Legend(at = c(-4, 0, max(4)), col_fun = color_scale, title = "centered.scaled")


lgd_list <- packLegend(ht_legend, lgd_species, lgd_subregion, direction = "vertical")


#png(here(plot_dir, "Heatmap_excitatory_celltypes.png"), width=10, height=6, units="in", res=300)
draw(ht, annotation_legend_list = lgd_list)
#dev.off()



