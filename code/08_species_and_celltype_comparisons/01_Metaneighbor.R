library("ggspavis")
library("scater")
library("pheatmap")
library("spatialLIBD")
library("patchwork")
library("scran")
library("Seurat")
library("SingleCellExperiment")
library(pheatmap)
library(MetaNeighbor)
library(ggplot2)
library(here)
library(ComplexHeatmap)

# directories
plot_dir = here("plots", "08_species_comparisons")

sce.excit <- readRDS(here("processed-data","07_annotation_and_characterization", "sce_excit_final_subclusters_annotated.rds"))
sce.inhib <- readRDS(here("processed-data","07_annotation_and_characterization", "sce_inhib_final_subclusters_annotated.rds"))


# ====== Copy necessary function from Metaneighbor =======

order_sym_matrix <- function(M, na_value = 0) {
    M <- (M + t(M))/2
    M[is.na(M)] <- na_value
    result <- stats::as.dendrogram(
        stats::hclust(stats::as.dist(1-M), method = "average")
    )
    return(result)
}

# ====== Inhibitory cells ======
#. set up colors
inhib.colors <- pals::cols25()[1:18]

inhib_celltype_colors <- setNames(inhib.colors, unique(sce.inhib$fine_celltype))

# we need to drop the CeA clusters; PRKCD_DRD1, PRKCD_DRD2, and ZFHX3_SCN5A
sce.inhib <- sce.inhib[, !names(sce.inhib$fine_celltype) %in% c("PRKCD_DRD1", "PRKCD_DRD2", "ZFHX3_SCN5A")]
inhib_celltype_colors <- inhib_celltype_colors[!names(inhib_celltype_colors) %in% c("PRKCD_DRD1", "PRKCD_DRD2", "ZFHX3_SCN5A")]

# reset levels
sce.inhib$fine_celltype <- factor(sce.inhib$fine_celltype)

# === Run US Metaneighbor ===

var_genes = variableGenes(dat = sce.inhib, exp_labels = sce.inhib$species)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = sce.inhib,
                             study_id = sce.inhib$species,
                             cell_type = sce.inhib$fine_celltype,
                             fast_version = TRUE,
                             one_vs_best=TRUE,
                             symmetric_output=FALSE)


#pdf(here(plot_dir,"meta_clusters.pdf"))
#plotMetaClusters(mclusters, celltype_NV)
#dev.off()

# === Heatmap Generation ===

# set color scale
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(100))
breaks = seq(0, 1, length=101)
ordering <- order_sym_matrix(celltype_NV)

# --- Heatmap Annotations ---
# get human, baboon, and macaque labels from the row names of celltype_NV
species_anno <- gsub("^(.*?)\\|.*$", "\\1", rownames(celltype_NV))

# get the 3rd rows celltype (CARTPT, PPP1R1B, SST) after | from the row names of celltype_NV
celltype_anno <- gsub("^.*?\\|(.*?)$", "\\1", rownames(celltype_NV))

# get distinct celltype annotation colors
# n <- length(celltype_anno)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]

# make a named vector from celltype_anno and col_vector
#names(col_vector) <- celltype_anno

# make top (column) annotations
row_ha <- rowAnnotation(
    species = species_anno,
    celltype = celltype_anno,
    col = list(
        species = c("baboon" ="#fe9380", "human" = "#b0d5f5", "macaque" = "#ffc84d"),
        celltype = inhib_celltype_colors
    ),
    show_annotation_name = FALSE,
    show_legend = c(TRUE, FALSE)
)

# make row annotations 
top_ha <- HeatmapAnnotation(
    celltype = celltype_anno,
    species = species_anno,
    col = list(
        species = c("baboon" ="#fe9380", "human" = "#b0d5f5", "macaque" = "#ffc84d"),
        celltype = inhib_celltype_colors
    ),
    show_annotation_name = FALSE,
    show_legend = c(TRUE, FALSE)
)

# get row order from heatmap to rename labels
hm <- ComplexHeatmap::Heatmap(celltype_NV,
                              cluster_rows = ordering,
                              cluster_columns = ordering,
                              show_column_dend = FALSE, 
                              show_row_dend = FALSE,
                              right_annotation=row_ha,
                              show_column_names=FALSE,
                              row_names_gp = grid::gpar(fontsize = 8),
                              name="AUROC",
                              show_row_names=FALSE
)
row_name_order <- row_order(hm)

# reorder celltype_anno based on row name order
celltype_labels <- celltype_anno[row_name_order]

# new row labels (only 1 for each celltype in the middle row)
new_row_labels <- rep("", length(celltype_labels))
for (i in seq(2, length(celltype_labels), 3)) {
    new_row_labels[i] <- celltype_labels[i]
}

# revert back to original row name order
new_row_labels <- new_row_labels[order(row_name_order)]

# plot complex heatmap
png(here(plot_dir, "Metaneighbor_heatmap_inhibitory_Samples.png"), width=10, height=6.5, units="in", res=300)
ComplexHeatmap::Heatmap(celltype_NV,
                        cluster_rows = ordering,
                        cluster_columns = ordering,
                        show_column_dend = FALSE, 
                        show_row_dend = FALSE,
                        right_annotation=row_ha,
                        top_annotation=top_ha,
                        show_column_names=FALSE,
                        row_names_gp = grid::gpar(fontsize = 8),
                        name="AUROC",
                        row_labels = new_row_labels,
                        column_title = "Inhibitory Celltypes",
                        column_title_side="top",
                        row_title=""
)
dev.off()



# ==== Identify top hits ===

celltype_NV_long <- as.data.frame(as.table(as.matrix(celltype_NV)))

# remove NAs
celltype_NV_long <- celltype_NV_long[!is.na(celltype_NV_long$Freq),]
celltype_NV_long

# Remove the species part from the Var2 column
data <- celltype_NV_long %>%
    mutate(Celltype2 = sub(".*\\|", "", Var2)) %>%
    mutate(Celltype1 = sub(".*\\|", "", Var1)) %>%
    select(Freq, Celltype1, Celltype2) %>%
    filter(Celltype1==Celltype2)

head(data)
# Freq    Celltype1    Celltype2
# 1 1.0000000 CARTPT_CDH23 CARTPT_CDH23
# 2 1.0000000 CARTPT_CDH23 CARTPT_CDH23
# 3 0.9988851 CARTPT_CDH23 CARTPT_CDH23
# 4 0.9954877     CCK_CNR1     CCK_CNR1
# 5 0.9874280     CCK_CNR1     CCK_CNR1
# 6 0.9815084     CCK_CNR1     CCK_CNR1


# Calculate mean AUROC for each cell type and sort them
sorted_celltypes <- data %>%
    group_by(Celltype1) %>%
    summarize(Mean_AUROC = mean(Freq)) %>%
    arrange(desc(Mean_AUROC)) %>%
    pull(Celltype1)

# Convert Celltype1 to a factor with levels ordered by Mean_AUROC
data$Celltype1 <- factor(data$Celltype1, levels = sorted_celltypes)

# Plot the boxplots using ggplot2
p<-ggplot(data, aes(x = Celltype1, y = Freq, fill = Celltype1)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cross-species celltype classification accuracy", y = "AUROC", x = "Celltype") +
    geom_hline(yintercept = 0.5, linetype = "dashed", size=1.5) +
    annotate("text", x = 1, y = 0.53, label = "Random chance", hjust = 0, size = 7) +
    # icnrease font sizes
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(size = 28)) +
    scale_fill_manual(values = inhib_celltype_colors) +
    theme(legend.position = "none") +
    xlab(NULL)

png(here(plot_dir, "Boxplots_AUROC_inhibitory.png"), width=10, height=6, units="in", res=300)
print(p)
dev.off()
    


# =========== Excitatory cells =============

# === Run US Metaneighbor ===
var_genes = variableGenes(dat = sce.excit, exp_labels = sce.excit$species)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = sce.excit,
                             study_id = sce.excit$species,
                             cell_type = sce.excit$fine_celltype,
                             fast_version = TRUE,
                             one_vs_best=TRUE,
                             symmetric_output=FALSE)


ordering <- order_sym_matrix(celltype_NV)
# === Heatmap Generation ===

excit.colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#FFA500", 
                  "#800080", "#008000", "#000080", "#FFC0CB", "#D2B48C")

excit_celltype_colors <- setNames(excit.colors, unique(sce.excit$fine_celltype))


# --- Heatmap Annotations ---
# get human, baboon, and macaque labels from the row names of celltype_NV
species_anno <- gsub("^(.*?)\\|.*$", "\\1", rownames(celltype_NV))

# get the 3rd rows celltype (CARTPT, PPP1R1B, SST) after | from the row names of celltype_NV
celltype_anno <- gsub("^.*?\\|(.*?)$", "\\1", rownames(celltype_NV))


# make top (column) annotations
row_ha <- rowAnnotation(
    species = species_anno,
    celltype = celltype_anno,
    col = list(
        species = c(human = "#AE93BEFF", baboon = "#B4DAE5FF", macaque = "#F0D77BFF"),
        celltype = excit_celltype_colors 
    ),
    show_annotation_name = FALSE,
    show_legend = c(TRUE, FALSE)
)

# make row annotations 
top_ha <- HeatmapAnnotation(
    celltype = celltype_anno,
    species = species_anno,
    col = list(
        species = c(human = "#AE93BEFF", baboon = "#B4DAE5FF", macaque = "#F0D77BFF"),
        celltype = excit_celltype_colors 
    ),
    show_annotation_name = FALSE,
    show_legend = c(TRUE, FALSE)
)

# get row order from heatmap to rename labels
hm <- ComplexHeatmap::Heatmap(celltype_NV,
                              cluster_rows = ordering,
                              cluster_columns = ordering,
                              show_column_dend = FALSE, 
                              show_row_dend = FALSE,
                              show_column_names=FALSE,
                              row_names_gp = grid::gpar(fontsize = 8),
                              name="AUROC",
                              show_row_names=FALSE
)
row_name_order <- row_order(hm)

# reorder celltype_anno based on row name order
celltype_labels <- celltype_anno[row_name_order]

# new row labels (only 1 for each celltype in the middle row)
new_row_labels <- rep("", length(celltype_labels))
for (i in seq(2, length(celltype_labels), 3)) {
    new_row_labels[i] <- celltype_labels[i]
}

# revert back to original row name order
new_row_labels <- new_row_labels[order(row_name_order)]

# plot complex heatmap
png(here(plot_dir, "Metaneighbor_heatmap_excitatory.png"), width=10, height=6.5, units="in", res=300)
ComplexHeatmap::Heatmap(celltype_NV,
                        cluster_rows = ordering,
                        cluster_columns = ordering,
                        show_column_dend = FALSE, 
                        show_row_dend = FALSE,
                        right_annotation=row_ha,
                        top_annotation=top_ha,
                        show_column_names=FALSE,
                        row_names_gp = grid::gpar(fontsize = 8),
                        name="AUROC",
                        row_labels = new_row_labels,
                        column_title = "Excitatory Celltypes"
)
dev.off()



# ==== Identify top hits ===

celltype_NV_long <- as.data.frame(as.table(as.matrix(celltype_NV)))

# remove NAs
celltype_NV_long <- celltype_NV_long[!is.na(celltype_NV_long$Freq),]
celltype_NV_long

# Remove the species part from the Var2 column
data <- celltype_NV_long %>%
    mutate(Celltype2 = sub(".*\\|", "", Var2)) %>%
    mutate(Celltype1 = sub(".*\\|", "", Var1)) %>%
    select(Freq, Celltype1, Celltype2) %>%
    filter(Celltype1==Celltype2)

head(data)
# Freq    Celltype1    Celltype2
# 1 1.0000000 CARTPT_CDH23 CARTPT_CDH23
# 2 1.0000000 CARTPT_CDH23 CARTPT_CDH23
# 3 0.9988851 CARTPT_CDH23 CARTPT_CDH23
# 4 0.9954877     CCK_CNR1     CCK_CNR1
# 5 0.9874280     CCK_CNR1     CCK_CNR1
# 6 0.9815084     CCK_CNR1     CCK_CNR1


# Calculate mean AUROC for each cell type and sort them
sorted_celltypes <- data %>%
    group_by(Celltype1) %>%
    summarize(Mean_AUROC = mean(Freq)) %>%
    arrange(desc(Mean_AUROC)) %>%
    pull(Celltype1)

# Convert Celltype1 to a factor with levels ordered by Mean_AUROC
data$Celltype1 <- factor(data$Celltype1, levels = sorted_celltypes)

# Plot the boxplots using ggplot2
p<-ggplot(data, aes(x = Celltype1, y = Freq, fill = Celltype1)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cross-species celltype classification accuracy", y = "AUROC", x = "Celltype") +
    geom_hline(yintercept = 0.5, linetype = "dashed", size=1.5) +
    annotate("text", x = 1, y = 0.53, label = "Random chance", hjust = 0, size = 7) +
    # icnrease font sizes
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(size = 28)) +
    scale_fill_manual(values = excit_celltype_colors) +
    theme(legend.position = "none") +
    xlab(NULL)

png(here(plot_dir, "Boxplots_AUROC_excitatory.png"), width=10, height=6, units="in", res=300)
print(p)
dev.off()



