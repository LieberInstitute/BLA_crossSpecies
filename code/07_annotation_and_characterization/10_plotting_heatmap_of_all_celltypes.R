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
library(ggh4x)

## save directories
plot_dir = here("plots", "07_annotation_and_characterization","01_broad_annotations")
#processed_dir = here("processed-data","07_annotation")
processed_dir <- here("processed-data")

# load sce
# sce.excit <- readRDS(here("processed-data", "07_annotation_and_characterization","sce_excit_final_subclusters_annotated.rds"))
# sce.inhib <- readRDS(here("processed-data","07_annotation_and_characterization", "sce_inhib_final_subclusters_annotated.rds"))
# sce.other <- readRDS(here("processed-data","07_annotation_and_characterization", "sce_other_final_celltypes.rds"))

# # combine SCE objects
# sce.inhib$subcluster_idents <- NULL
# sce.inhib$integrated_snn_res.0.4 <- NULL
# sce.excit$integrated_snn_res.0.4 <- NULL

# sce.neurons <- cbind(sce.excit, sce.inhib)
# sce <- cbind(sce.neurons, sce.other)

# rm(sce.excit, sce.inhib, sce.other, sce.neurons)
# gc()

# save sce
# saveRDS(sce, here("processed-data","07_annotation_and_characterization","sce_FINAL_all_celltypes.rds"))

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


# ========= Heatmap ==========

genes <- c(
    "SNAP25", 
    "SYT1",
    "SLC17A7",
    "SLC17A6",
    "GAD1",
    "GAD2"
    #"SIM1",
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
# species_colors <- c("baboon" ="#fe9380", "human" = "#b0d5f5", "macaque" = "#ffc84d")
# subregion_colors <- c( "Accessory Basal" ="#8A2BE2", "Basal" = "#FF4500","Central Nucleus" = "#FFD700", "Lateral" = "#00BFFF")




# cell_types <- sce$fine_celltype
# species <- sce$species

# bar_width <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

# species_proportions <- table(cell_types, species) / rowSums(table(cell_types, species))
# ha_species <- rowAnnotation(Species = anno_barplot(species_proportions, gp = gpar(fill = species_colors),
#                                                    bar_width=.8
# ), 
# border=FALSE,
# width=unit(2,"cm")
# )



# == Creating stacked bar plot subregion annotations == 
# sce.macaque <- sce[,which(colData(sce)$species == "macaque")]
# cell_types <- sce.macaque$fine_celltype
# subregions <- sce.macaque$Subregion  

# bar_width <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

# subregion_proportions <- table(cell_types, subregions) / rowSums(table(cell_types, subregions))
# ha_subregion <- rowAnnotation(Subregion = anno_barplot(subregion_proportions, gp = gpar(fill = subregion_colors),
#                                                        bar_width=.8
# ), 
# border=FALSE,
# width=unit(2,"cm"),
# annotation_legend_param = list(title = "Subregion Proportions")
# )



# #Heatmap Scale
color_scale <- circlize::colorRamp2(c(-4, 0, max(scaled_mat)), c("blue", "#EEEEEE", "red"))
row_split <- factor(c(rep("Excitatory", 12), rep("Inhibitory", 18), rep("Non-neuronal", 6)))
# Manually define the row order based on the heatmap you provided
row_order <- c(
  "ADARB2_TRPS1", "ZBTB20_SLC4A4", "GULP1_TRHDE", "PEX5L_MYRIP", "ESR1_ADRA1A", 
  "MEIS2_COL25A1", "GRIK3_TNS3", "ST18_ABCA8", "SATB2_MPPED1", "RXFP1_KIAA1217", 
  "SLC17A8_ST8SIA2",  "MEIS1_PARD3B","CARTPT_CDH23", "ST18_IL1RAPL2", "LAMP5_NTNG1", 
  "LAMP5_KIT", "PVALB_MYO5B", "SST_PRKCQ", "CCK_CNR1", "VIP_PLPP4", 
  "VIP_ADRA1B", "THSD7B_CALB2", "ZFHX3_SCN5A", "LHX8_ANGPT1", "TSHZ1.2", 
  "TSHZ1.1", "PRKCD_DRD2", "PRKCD_DRD1", "PVALB_UNC5B", "SST_NOS1", 
  "Endothelial", "Astrocyte", "Ependymal", "Oligodendrocyte", "OPC",  "Microglia"
)

# Match the desired row names to the actual row names in your matrix
row_order_indices <- match(row_order, rownames(scaled_mat))

scaled_mat <- scaled_mat[row_order, ]

# == Creating celltype label annotations ===
# Create a named vector mapping cell types to colors
# 
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

# save celltype colors
saveRDS(celltype_colors, here("processed-data", "07_annotation_and_characterization", "celltype_colors.rds"))

anno_df <- data.frame(Cell_Type = out$fine_celltype)
anno_df <- anno_df[row_order_indices, , drop = FALSE]

ha_celltypes <- rowAnnotation(Cell_Type = anno_df$Cell_Type, 
                              col = list(Cell_Type = celltype_colors),
                              annotation_legend_param = list(Cell_Type = list(title = "Cell Types")),
                              show_legend=FALSE
)


ht <- Heatmap(scaled_mat,
              row_split = row_split, # Specify the custom row splits
              cluster_row_slices = FALSE,
              cluster_columns=FALSE,
              col=color_scale,
              rect_gp = gpar(col = "black", lwd = .75),
              left_annotation=ha_celltypes,
              show_heatmap_legend = TRUE,
              row_names_side="left",
              show_row_dend = FALSE,
              row_order=row_order,
              row_title_gp = gpar(fontface = "bold"),  # Bold the row titles
              row_gap = unit(2, "mm"),
                           heatmap_legend_param = list(title = "Lognorm Counts\nCentered, Scaled", 
                                           title_position = "leftcenter",
                                           direction = "horizontal", 
                                           legend_height = unit(2, "cm"))  # Adjust height as needed

) 
  #  ha_species +
  #  ha_subregion

#ht_legend <- Legend(at = c(-4, 0, max(4)), col_fun = color_scale, title = "centered.scaled", direction="horizontal")

png(here(plot_dir, "Heatmap_excitatory_celltypes_NEW.png"), width=4, height=9, units="in", res=300)
draw(ht, heatmap_legend_side = "top")
dev.off()




# ==============================================
#  Proportion stacked bar plots for Figure 1 
# ==============================================

library(dplyr)
# ======= Macaque samples ========

#subset to different species
sce_macaque <- sce[,which(colData(sce)$species == "macaque")]

# drop non-neuronal
sce_macaque <- sce_macaque[,which(colData(sce_macaque)$broad_celltype != "Non-neuronal")]

# Create the table
celltype_table <- table(sce_macaque$fine_celltype, sce_macaque$Sample, sce_macaque$Subregion)

# Convert the table to a data frame
df <- as.data.frame(celltype_table)

# Rename columns for clarity
colnames(df) <- c("CellType", "Sample", "Subregion", "Count")

# Calculate proportions within each Sample and Subregion
df <- df %>%
  group_by(Sample, Subregion) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  filter(!is.na(Proportion))

# change Sample names to 1 - max
df$Sample <- as.character(as.numeric(as.factor(df$Sample)))

# change Accessory Basal to Acc. Basal
levels(df$Subregion)[levels(df$Subregion) == "Accessory Basal"] <- "acBasal"
levels(df$Subregion)[levels(df$Subregion) == "Central Nucleus"] <- "Central"


# ==== Define custom colors for cell types ====

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
# sce.other <- sce[,which(colData(sce)$broad_celltype == "Other")]
# greys <- gray.colors(6)
# other_colors <- setNames(greys, unique(sce.other$fine_celltype))

# combine all colors
celltype_colors <- c(inhibitory_colors, excitatory_colors)



# ==== Plot =====
 # Create the stacked bar plot with facet wrap for subregion
png(here(plot_dir, "AllNuclei_stackedbars_macaque_sample_by_celltype_subregion.png"), width = 10, height = 5, units = "in", res = 300)
p1 <- ggplot(df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "Proportion", fill = "Cell Type") +
  theme(axis.text.x = element_text(size=12),
        strip.text = element_text(size = 18),
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=18)) +  # Customize facet labels
  scale_fill_manual(values = celltype_colors) +
  facet_grid(~Subregion, scales="free", space="free") +  # Adjust the number of columns in the facet wrap
    theme(legend.position = "none")
p1
dev.off()




# ======= Macaque samples ========


sce.neun <- sce[,which(colData(sce)$broad_celltype != "Non-neuronal")]

# Replace NA values in Subregion with "Whole BLA"
sce.neun$Subregion <- as.character(sce.neun$Subregion)
sce.neun$Subregion[is.na(sce.neun$Subregion)] <- "Whole BLA"
sce.neun$Subregion <- factor(sce.neun$Subregion)


# Create the table
celltype_table <- table(sce.neun$fine_celltype, sce.neun$Sample, sce.neun$Subregion, sce.neun$species)

# Convert the table to a data frame
df <- as.data.frame(celltype_table)

# Rename columns for clarity
colnames(df) <- c("CellType", "Sample", "Subregion", "Species", "Count")

# Convert Subregion to character before modifying it
# df$Subregion <- as.character(df$Subregion)
# df$Subregion[df$Species == "human"] <- "Basal"
# df$Subregion <- factor(df$Subregion)

# Calculate proportions within each Sample and Subregion
df <- df %>%
  group_by(Species, Sample, Subregion) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  filter(!is.na(Proportion))



# change Sample names to 1 - max
df$Sample <- as.character(as.numeric(as.factor(df$Sample)))

# change Accessory Basal to Acc. Basal
levels(df$Subregion)[levels(df$Subregion) == "Accessory Basal"] <- "acBasal"
levels(df$Subregion)[levels(df$Subregion) == "Central Nucleus"] <- "Central"


# capatilize species names
df$Species <- gsub("human", "Human", df$Species)
df$Species <- gsub("macaque", "Macaque", df$Species)
df$Species <- gsub("baboon", "Baboon", df$Species)


# ==== Define custom colors for cell types ====

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
# sce.other <- sce[,which(colData(sce)$broad_celltype == "Other")]
# greys <- gray.colors(6)
# other_colors <- setNames(greys, unique(sce.other$fine_celltype))

# combine all colors
celltype_colors <- c(inhibitory_colors, excitatory_colors)



# ==== Plot =====
# Reorder Species levels (adjust the order as per your requirement)
df$Species <- factor(df$Species, levels = c("Human", "Baboon", "Macaque"))  # Replace with your actual species names in desired order

# reorder Subregion levels to Whole BLA, Lateral, Basal, Accessory Basal, Central Nucleus
df$Subregion <- factor(df$Subregion, levels = c("Whole BLA", "Lateral", "Basal", "acBasal", "Central"))

# Create the stacked bar plot with facet wrap for subregion
png(here(plot_dir, "AllNuclei_stackedbars_celltype_by_subregion.png"), width = 24, height = 7.5, units = "in", res = 300)
p1 <- ggplot(df, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Proportion", fill = "Cell Type") +
  scale_fill_manual(values = celltype_colors) +
  facet_nested(.~Species+Subregion, scales = "free", space = "free", drop=TRUE) +  # Adjust the number of columns in the facet wrap
    theme_gray() +
  theme(axis.text.x = element_text(size=22),
        strip.text = element_text(size = 28), # Customize facet labels 
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=28),  
        axis.title.x=element_text(size=28),   
        legend.position = "none") 

p1
dev.off()






# =========== New UMAPs for Figure 1 ===========
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

region_colors <- c("Basal" = "#E41A1C", "Lateral" = "#377EB8", "Central Nucleus" = "#4DAF4A", "Accessory Basal" = "#984EA3", "NA" = "black")

# reorder again to best show differences across subregions
sce$species <- forcats::fct_relevel(sce$species, "human", "baboon", "macaque")
sce <- sce[,order(colData(sce)$species)]

png(here(plot_dir, "UMAPs_Figure1.png"), width=15, height=5, units="in", res=300)

p1 <- plotReducedDim(sce, dimred = "UMAP_new", colour_by = "Subregion", point_size=0.2) +
    scale_color_manual(values = region_colors) +
    theme_void() +
    ggtitle("Subregion")  +
    theme(axis.line = element_blank(), # Remove the default axis lines
    axis.text = element_blank(), # Remove axis text
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(size = 18),
    legend.position = c(0.82, 0.85),
    legend.text = element_text(size = 12), # Increase legend text size
    legend.title = element_blank()
   ) +
   guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))


p2 <- plotReducedDim(sce, dimred = "UMAP_new", colour_by = "broad_celltype", point_size=0.2) +
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


p3 <- plotReducedDim(sce, dimred = "UMAP_new", colour_by = "fine_celltype", point_size=0.2) +
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

p1+p2+p3
dev.off()

