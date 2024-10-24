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
library(dplyr)


## save directories
plot_dir = here("plots", "07_annotation_and_characterization","05_inhib_annotations")
processed_dir <- here("processed-data")

# load sce
sce <- readRDS(here("processed-data","07_annotation_and_characterization", "sce_inhib_final_subclusters_annotated.rds"))
sce

# ======= UMAPs ========


#subset to different species
sce_macaque <- sce[,which(colData(sce)$species == "macaque")]
sce_baboon <- sce[,which(colData(sce)$species == "baboon")]
sce_human <- sce[,which(colData(sce)$species == "human")]

dim(sce)
# [1] 13842 29606

dim(sce_macaque)
# [1] 13842 18724

dim(sce_baboon)
# [1] 13842  6464

dim(sce_human)
# [1] 13842  4418

# Define the number of colors you want
# mycolors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#FFA500", 
#             "#800080", "#008000", "#000080", "#FFC0CB", "#D2B48C", "#808000", "#000000", 
#             "#808080", "#FF1493", "#1E90FF", "#8B4513")




# ========== UMAP for top of figure ==========

mycolors <- pals::cols25()[1:18]
celltype_colors <- setNames(mycolors, unique(sce$fine_celltype))
celltype_colors

# == calculating centroid for geom_labels ==
# Extract t-SNE coordinates
umap_coords <- reducedDim(sce, "UMAP")

# Create a data frame for plotting
plot_data <- data.frame(umap_coords, Celltype = sce$fine_celltype)

# Calculate centroids for each group
centroid_data <- plot_data %>%
  group_by(Celltype) %>%
  summarize(CentroidX = mean(UMAP_1), CentroidY = mean(UMAP_2), .groups = 'keep')

# Match colors to annotations in centroid_data
centroid_data$Color <- celltype_colors[centroid_data$Celltype]


png(here(plot_dir, "UMAP_inhib_annotated_celltype_colors.png"), width=7, height=7, units="in", res=300)
p1 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "fine_celltype", point_size=0.75) +
    scale_color_manual(values = celltype_colors) +
  ggrepel::geom_label_repel(data = centroid_data, aes(x = CentroidX, y = CentroidY, label = Celltype),
               #color = centroid_data$Color,  
               fontface = "bold.italic", size=5) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p1
dev.off()
  

# ========= Plotting UMAPs with subregions =========

unique(sce_baboon$Subregion)
# [1] "Basal"   "Lateral"

unique(sce_macaque$Subregion)
# [1] "Lateral"         "Basal"           "Central Nucleus" "Accessory Basal"

# set colors so that they are consistent across species
subregion_colors <- c( "Accessory Basal" ="#8A2BE2", "Basal" = "#FF4500","Central Nucleus" = "#FFD700", "Lateral" = "#00BFFF")


# plot UMAP with idents
png(here(plot_dir, "UMAP_inhib_across_species_subregions.png"), width=12, height=4, units="in", res=300)

p1 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = subregion_colors) +
    theme_void() +
    ggtitle("All nuclei (n = 29,606)")  +
    theme(legend.position="none")

p2 <- plotReducedDim(sce_macaque, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = subregion_colors) +
    theme_void() +
    ggtitle("Macaque (n = 18,724)")  +
    theme(legend.position="none")

p3 <- plotReducedDim(sce_baboon, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = subregion_colors) +
    theme_void() +
    ggtitle("Baboon (n = 6,464)") +
    theme(legend.position="none")

p4 <- plotReducedDim(sce_human, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    scale_color_manual(values = subregion_colors) +
    theme_void() +
    ggtitle("Human (n = 4,418)")  +
    theme(legend.position="none")


ggarrange(p1, p2, p3, p4, ncol=4, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()



# ======================================
#  Heatmap 
# ======================================

genes <- c("CARTPT","CCK", "LAMP5","PRKCD", "DRD1","DRD2","PENK", "PVALB", "NPY","SST","NOS1",
           "CALB1","CALB2", "FOXP2", "GAL", "VIP","CRH","NTS")

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
species_colors <- c("baboon" ="#1b4543", "human" = "#f0be6f", "macaque" = "#b3d0c6")



cell_types <- sce$fine_celltype
species <- sce$species

bar_width <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

# Calculate the total number of cells in each subregion
total_cells_per_species <- table(species)

# Calculate the desired proportion 
desired_proportion <- 1/3

# Calculate the scaling factor for each subregion
scaling_factors <- desired_proportion / (total_cells_per_species / sum(total_cells_per_species))

# Create a contingency table of cell types vs. species
cell_type_species_table <- table(cell_types, species)

# Adjust the cell counts by scaling factor
adjusted_cell_counts <- sweep(cell_type_species_table, 2, scaling_factors, "*")

# Calculate the adjusted proportions
adjusted_proportions <- adjusted_cell_counts / rowSums(adjusted_cell_counts)

#species_proportions <- table(cell_types, species) / rowSums(table(cell_types, species))
ha_species <- rowAnnotation(Species = anno_barplot(adjusted_proportions, gp = gpar(fill = species_colors),
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


# Calculate the total number of cells in each subregion
total_cells_per_subregion <- table(subregions)

# Calculate the desired proportion (25% for each of the 4 subregions)
desired_proportion <- 0.25

# Calculate the scaling factor for each subregion
scaling_factors <- desired_proportion / (total_cells_per_subregion / sum(total_cells_per_subregion))

# Create a contingency table of cell types vs. subregions
cell_type_subregion_table <- table(cell_types, subregions)

# Adjust the cell counts by scaling factor
adjusted_cell_counts <- sweep(cell_type_subregion_table, 2, scaling_factors, "*")

# Calculate the adjusted proportions
adjusted_proportions <- adjusted_cell_counts / rowSums(adjusted_cell_counts)

#subregion_proportions <- table(cell_types, subregions) / rowSums(table(cell_types, subregions))
ha_subregion <- rowAnnotation(Subregion = anno_barplot(adjusted_proportions, gp = gpar(fill = subregion_colors),
                                                       bar_width=.8
                                                       ), 
                              border=FALSE,
                              width=unit(2,"cm"),
                              annotation_legend_param = list(title = "Subregion Proportions")
                              )
                              

# == Creating celltype label annotations ===
# Create a named vector mapping cell types to colors

# colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#FFA500", 
#             "#800080", "#008000", "#000080", "#FFC0CB", "#D2B48C", "#808000", "#000000", 
#             "#808080", "#FF1493", "#1E90FF", "#8B4513")


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
                        legend_gp = gpar(fill = subregion_colors, fontsize=16), pch = 16)



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


png(here(plot_dir, "Heatmap_inhibitory_celltypes.png"), width=10, height=6, units="in", res=300)
draw(ht, annotation_legend_list = lgd_list)
dev.off()





# ============ Porportion of cell-tyoes per putative region ========


# Filter the data for macaque species
sce.macaque <- sce[, which(colData(sce)$species == "macaque")]
cell_types <- sce.macaque$fine_celltype
subregions <- sce.macaque$Subregion  

bar_width <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

# Create contingency table
contingency_table <- table(cell_types, subregions)

# Calculate the proportion of each cell type per subregion
subregion_proportions <- prop.table(contingency_table, margin = 2)

# Convert the data to long format
data_long <- reshape2::melt(subregion_proportions, variable.name = "Subregion", value.name = "Proportion")
colnames(data_long) <- c("CellType", "Subregion", "Proportion")

# Set the levels of Subregion factor in the desired order
data_long$Subregion <- factor(data_long$Subregion, levels = c("Lateral", "Basal","Accessory Basal", "Central Nucleus"))

# Create the stacked bar plot
p1 <- ggplot(data_long, aes(x = Subregion, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    labs(title = "Proportion of inhibitory cell-types in Macaque",
         y = "Proportion",
         x=NULL) +
    scale_fill_manual(name = "Cell Type", values = celltype_colors) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #theme(legend.position = "none") +
    coord_flip() +
    # increase font sizes
    theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(size = 20, hjust = 0.5))

png(here(plot_dir, "Proportion_inhibitory_celltypes_per_subregion.png"), width=10, height=4.5, units="in", res=300)
p1
dev.off()





# ============ Tukey test for significant differences in cell-type proportions =========

library(rstatix)
library(dplyr)
library(ggpubr)

# Filter the data for macaque species
sce.macaque <- sce[, which(colData(sce)$species == "macaque")]

# Extract cell types, subregions, and sample information
cell_types <- sce.macaque$fine_celltype
subregions <- sce.macaque$Subregion  
samples <- sce.macaque$Sample

# Create a data frame with counts per sample
cell_data <- data.frame(
  CellType = cell_types,
  Subregion = subregions,
  Sample = samples
)

# Calculate the proportion of each cell type per sample
cell_proportions <- cell_data %>%
  group_by(Sample, Subregion, CellType) %>%
  summarize(Count = n(), .groups = "drop") %>%  # Add .groups = "drop"
  group_by(Sample, Subregion) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()


# Perform post-hoc Tukey test for each cell type and store results
tukey_results <- cell_proportions %>%
  group_by(CellType) %>%
  tukey_hsd(Proportion ~ Subregion) %>%
  add_significance() %>%
  filter(p.adj.signif != "ns") %>%
  add_xy_position(x="Subregion")  # Filter to keep only significant results

# Modify the y.position for plotting (for manual adjustment of p-value locations)
tukey_results <- tukey_results 

# Now plot using ggplot2 and manually add the significant p-values
p3 <- ggplot(cell_proportions, aes(x = Subregion, y = Proportion)) +
  geom_boxplot(aes(fill = Subregion)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  #theme_minimal() +
  labs(title = "Proportion of Cell Types per Sample by Subregion in Macaques", 
       y = "Proportion", 
       x = "Subregion") +
  scale_fill_manual(name = "Subregion", values = subregion_colors) +
  facet_wrap(~ CellType, scales = "free") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.text = element_text(size = 14),
        # Adjust legend size
        legend.text = element_text(size = 14),  # Increase legend text size
        legend.title = element_text(size = 16),  # Increase legend title size
        legend.key.size = unit(1.5, "lines")) +  # Increase legend key size) +
  stat_pvalue_manual(tukey_results, label.size=6)

# Save the plot with the significant p-values only
png(here(plot_dir, "Faceted_Boxplot_CellType_Subregion_Proportion_SignificantOnly.png"), width = 13, height = 13, units = "in", res = 300)
p3
dev.off()






library(rstatix)
library(dplyr)
library(ggpubr)

# Filter the data for macaque species
sce.macaque <- sce[, which(colData(sce)$species == "macaque")]

# Extract cell types, subregions, and sample information
cell_types <- sce.macaque$fine_celltype
subregions <- sce.macaque$Subregion  
samples <- sce.macaque$Sample

# Create a data frame with counts per sample
cell_data <- data.frame(
  CellType = cell_types,
  Subregion = subregions,
  Sample = samples
)

# List of samples to relabel and highlight
samples_to_relabel <- c("VC_snRNAseq-7_LateralDorsalAP1AP2_-8", "VC_snRNAseq_8_Animal3_LD")

# Modify the Subregion for those specific samples and keep the rest unchanged
cell_data <- cell_data %>%
  mutate(Subregion = as.character(Subregion)) %>%
  mutate(Subregion = ifelse(Sample %in% samples_to_relabel, 
                            "Central Nucleus", Subregion))  # Only change for these specific samples

# Recalculate the cell proportions after modifying the subregions
cell_proportions <- cell_data %>%
  group_by(Sample, Subregion, CellType) %>%
  summarize(Count = n(), .groups = "drop") %>%
  group_by(Sample, Subregion) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Add a new column to highlight swapped samples (set as "Swapped" or "Normal")
cell_proportions <- cell_proportions %>%
  mutate(Relabel = ifelse(Sample %in% samples_to_relabel, "LA > CeA", "Clean"))

# Set the levels of Subregion factor in the desired order
cell_proportions$Subregion <- factor(cell_proportions$Subregion, 
                                     levels = c("Lateral", "Basal", "Accessory Basal", "Central Nucleus"))

# Perform post-hoc Tukey test for each cell type and store results
tukey_results <- cell_proportions %>%
  group_by(CellType) %>%
  tukey_hsd(Proportion ~ Subregion) %>%
  add_significance() %>%
  add_xy_position(x = "Subregion", step.increase = 0.1)  %>%
  filter(p.adj.signif != "ns") 


# Define colors for all subregions
subregion_colors <- c(
  "Lateral" = "#00BFFF", 
  "Basal" = "#FF4500", 
  "Accessory Basal" = "#8A2BE2", 
  "Central Nucleus" = "#FFD700",
  "Swapped" = "#0000FF"  # Blue for the swapped samples
)

# Now plot using ggplot2 and manually add the significant p-values
p3 <- ggplot(cell_proportions, aes(x = Subregion, y = Proportion)) +
  geom_boxplot(aes(fill=Subregion)) +
  geom_jitter(aes(color = Relabel), width = 0.2, alpha = 0.5) +  # Color the swapped samples differently
  labs(title = "Proportion of Cell Types per Sample by Subregion in Macaques", 
       y = "Proportion", 
       x = "Subregion") +
  scale_fill_manual(name = "Subregion", values = subregion_colors) +  # Explicitly set colors for each subregion
  scale_color_manual(values = c("LA > CeA" = "#00BFFF")) +  # Highlight swapped samples in blue
  facet_wrap(~ CellType, scales = "free") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),  # Increase legend text size
        legend.title = element_text(size = 16),  # Increase legend title size
        legend.key.size = unit(1.5, "lines")) +  # Increase legend key size
  stat_pvalue_manual(tukey_results, label.size = 6, , label = "p.adj.signif")

# Save the plot with the significant p-values only
png(here(plot_dir, "Faceted_Sample_Swapped_Boxplot_CellType_Subregion_Proportion_SignificantOnly.png"), width = 13, height = 13, units = "in", res = 300)
print(p3)
dev.off()  










# Define the subset of cell types you're interested in
selected_celltypes <- c("CCK_CNR1", "LAMP5_KIT", "PRKCD_DRD1", "PRKCD_DRD2", 
                        "PVALB_MYOSB", "SST_NOS1", "SST_PRKCQ", "TSHZ1.1", "VIP_ADRA1B")

# Filter cell_proportions and tukey_results for the selected cell types
cell_proportions_subset <- cell_proportions %>%
  filter(CellType %in% selected_celltypes)

tukey_results_subset <- tukey_results %>%
  filter(CellType %in% selected_celltypes)

# Now plot using ggplot2 and manually add the significant p-values
p_subset <- ggplot(cell_proportions_subset, aes(x = Subregion, y = Proportion)) +
  geom_boxplot(aes(fill = Subregion)) +
  geom_jitter(aes(color = Relabel), width = 0.2, alpha = 0.5) +  # Color the swapped samples differently
  labs(title = "Proportion of Cell Types across Subregion", 
       y = "Proportion", 
       x = "Subregion") +
  scale_fill_manual(name = "Subregion", values = subregion_colors) +  # Explicitly set colors for each subregion
  scale_color_manual(values = c("LA > CeA" = "#00BFFF")) +  # Highlight swapped samples in blue
  facet_wrap(~ CellType, nrow = 2, scales = "free") +  # Plot on two rows
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),  # Increase legend text size
        legend.title = element_text(size = 16),  # Increase legend title size
        legend.key.size = unit(1.5, "lines")) +  # Increase legend key size
  stat_pvalue_manual(tukey_results_subset, label.size = 6, label = "p.adj.signif")

# Save the plot with the significant p-values only
png(here(plot_dir, "Selected_Celltypes_Boxplot_TwoRows.png"), width = 10, height = 6, units = "in", res = 300)
print(p_subset)
dev.off()
