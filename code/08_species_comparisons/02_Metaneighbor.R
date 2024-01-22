library(dplyr)
library("Matrix")
library("here")
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
library(tidyverse)
library(ggplot2)
library(tidySingleCellExperiment)
library(ghibli)

# directories
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "08_species_comparisons")

sce.inhib <- readRDS(here(processed_dir, "sce.inhib.final.rds"))
sce.excit <- readRDS(here(processed_dir, "sce.excit.integrated.annotated.rds"))


# ====== Inhibitory cells ======

# === Run US Metaneighbor ===

var_genes = variableGenes(dat = sce.inhib, exp_labels = sce.inhib$species)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = sce.inhib,
                             study_id = sce.inhib$species,
                             cell_type = sce.inhib$fine_type,
                             fast_version = TRUE,
                             one_vs_best=TRUE,
                             symmetric_output=FALSE)


pdf(here(plot_dir,"meta_clusters.pdf"))
plotMetaClusters(mclusters, celltype_NV)
dev.off()

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
n <- length(celltype_anno)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]

# make a named vector from celltype_anno and col_vector
names(col_vector) <- celltype_anno

# make top (column) annotations
row_ha <- rowAnnotation(
    species = species_anno,
    celltype = celltype_anno,
    col = list(
        species = c(human = "#AE93BEFF", baboon = "#B4DAE5FF", macaque = "#F0D77BFF"),
        celltype = col_vector
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
        celltype = col_vector
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
pdf(here(plot_dir, "Metaneighbor_heatmap_inhibitory.pdf"), width=10, height=6.5)
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
                        column_title_side="bottom",
                        row_title=""
)
dev.off()

# =========== Excitatory cells =============

# === Run US Metaneighbor ===

var_genes = variableGenes(dat = sce.excit, exp_labels = sce.excit$species)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = sce.excit,
                             study_id = sce.excit$species,
                             cell_type = sce.excit$fine_type,
                             fast_version = TRUE,
                             one_vs_best=TRUE,
                             symmetric_output=FALSE)



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
n <- length(celltype_anno)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]

# make a named vector from celltype_anno and col_vector
names(col_vector) <- celltype_anno

# make top (column) annotations
row_ha <- rowAnnotation(
    species = species_anno,
    celltype = celltype_anno,
    col = list(
        species = c(human = "#AE93BEFF", baboon = "#B4DAE5FF", macaque = "#F0D77BFF"),
        celltype = col_vector
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
        celltype = col_vector
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
pdf(here(plot_dir, "Metaneighbor_heatmap_excitatory.pdf"), width=10, height=6.5)
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
# ========== Identify top hits ==========

top_hits = topHits(cell_NV = celltype_NV,
                   dat = sce.inhib,
                   study_id = sce.inhib$species,
                   cell_type = sce.inhib$fine_type,
                   threshold = 0.9)
top_hits
Study_ID|Celltype_1  Study_ID|Celltype_2 Mean_AUROC         Match_type
# 1      baboon|SST_NOS1     macaque|SST_NOS1       1.00 Reciprocal_top_hit
# 2        human|PVALB.3       baboon|PVALB.3       1.00 Reciprocal_top_hit
# 3        human|PVALB.3      macaque|PVALB.3       1.00          Above_0.9
# 4        human|PVALB.2       baboon|PVALB.2       1.00 Reciprocal_top_hit
# 5       baboon|PVALB.2      macaque|PVALB.2       1.00          Above_0.9
# 6     macaque|SST_NOS1       human|SST_NOS1       1.00          Above_0.9
# 7            human|CCK           baboon|CCK       1.00 Reciprocal_top_hit
# 8            human|CCK          macaque|CCK       1.00          Above_0.9
# 9       baboon|LAMP5.2      macaque|LAMP5.2       1.00 Reciprocal_top_hit
# 10  human|PPP1R1B_DRD1  baboon|PPP1R1B_DRD1       1.00 Reciprocal_top_hit
# 11          baboon|SST          macaque|SST       1.00 Reciprocal_top_hit
# 12     human|PENK_DRD2    macaque|PENK_DRD2       0.99 Reciprocal_top_hit
# 13      baboon|LAMP5.1      macaque|LAMP5.1       0.99 Reciprocal_top_hit
# 14     macaque|LAMP5.2        human|LAMP5.2       0.99          Above_0.9
# 15          baboon|SST            human|SST       0.99          Above_0.9
# 16       human|TSHZ1.1       baboon|TSHZ1.1       0.99 Reciprocal_top_hit
# 17     macaque|LAMP5.1        human|LAMP5.1       0.99          Above_0.9
# 18      baboon|TSHZ1.1      macaque|TSHZ1.1       0.99          Above_0.9
# 19     human|PENK_DRD2     baboon|PENK_DRD2       0.99          Above_0.9
# 20         human|VIP.2         baboon|VIP.2       0.99 Reciprocal_top_hit
# 21         human|VIP.2        macaque|VIP.2       0.99          Above_0.9
# 22      baboon|TSHZ1.2      macaque|TSHZ1.2       0.98 Reciprocal_top_hit
# 23       human|PVALB.1      macaque|PVALB.1       0.98 Reciprocal_top_hit
# 24       human|PVALB.1       baboon|PVALB.1       0.98          Above_0.9
# 25  human|PPP1R1B_DRD1 macaque|PPP1R1B_DRD1       0.98          Above_0.9
# 26      baboon|TSHZ1.2        human|TSHZ1.2       0.98          Above_0.9
# 27       baboon|CARTPT       macaque|CARTPT       0.97 Reciprocal_top_hit
# 28      macaque|CARTPT         human|CARTPT       0.96          Above_0.9
# 29          human|RELN          baboon|RELN       0.96 Reciprocal_top_hit
# 30 baboon|PPP1R1B_DRD1          human|ZFHX3       0.96          Above_0.9
# 31        baboon|VIP.1        macaque|VIP.1       0.96 Reciprocal_top_hit
# 32         baboon|RELN         macaque|RELN       0.95          Above_0.9
# 33        baboon|VIP.1          human|VIP.1       0.95          Above_0.9
# 34         human|ZFHX3        macaque|ZFHX3       0.95          Above_0.9
# 35         human|ZFHX3         baboon|ZFHX3       0.94          Above_0.9




pdf(here(plot_dir,"Metaneighbor_AUROC_interneurons_shuffled.pdf"))
AUROC_scores_shuffled = MetaNeighbor(dat = sce.inhib,
                                     experiment_labels =  as.numeric(factor(sce.inhib$species)),
                                     celltype_labels = celltype_matrix,
                                     genesets = var_genes,
                                     bplot = TRUE,
                                     fast_version = TRUE)
dev.off()











celltype_matrix = model.matrix(~sce.inhib$fine_type - 1)
colnames(celltype_matrix) = levels(as.factor(sce.inhib$fine_type))

# load in GO terms
#data(GOhuman)
GOhuman = readRDS("go_human.rds") 
length(GOhuman)
# [1] 22963

known_genes = rownames(sce.inhib)
go_sets = lapply(GOhuman, function(gene_set) { gene_set[gene_set %in% known_genes] }) 
min_size = 10
max_size = 100
go_set_size = sapply(go_sets, length)
go_sets = go_sets[go_set_size >= min_size & go_set_size <= max_size]
length(go_sets)
# [1] 6695

pdf(here(plot_dir,"Metaneighbor_AUROC_interneurons_extended.pdf"))
AUROC_scores = MetaNeighbor(dat = celltypes_of_interest ,
                            experiment_labels = as.numeric(factor(celltypes_of_interest$study_id)),
                            celltype_labels = celltype_matrix,
                            genesets = GOhuman,
                            bplot = TRUE,
                            fast_version = TRUE)
dev.off()

# plotting using ggplot
long_data <- data.frame(AUROC_scores) %>%
    rownames_to_column(var = "GO_term") %>%
    gather(Celltype, AUROC, -GO_term)
head(long_data) 
# GO_term   Celltype AUROC
# 1 GO:0000228 HTR3A.DRD2 0.915
# 2 GO:0000902 HTR3A.DRD2 0.996
# 3 GO:0000988 HTR3A.DRD2 0.964
# 4 GO:0003013 HTR3A.DRD2 0.991
# 5 GO:0003729 HTR3A.DRD2 0.848
# 6 GO:0003735 HTR3A.DRD2 0.657

ggplot(long_data, aes(x = Celltype, y = AUROC)) +
    geom_violin(aes(fill = Celltype), trim = FALSE) +  # Violin plot
    geom_jitter(width = 0.2, size = 2, alpha = 0.3) +  # Individual data points
    stat_summary(aes(group = interaction(Celltype)), fun = mean, geom = "crossbar", color = "inhibck", width=.5, position = position_dodge(width = 0.8)) +
    labs(title = "AUROC scores of interneurons using GO terms", y = "AUROC", x = "Celltype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
    coord_cartesian(ylim = c(0, 1))  # Set y-axis limits to 0 and 1

ggsave(here(plot_dir, "Metaneighbor_AUROC_interneuron_matched.pdf"), width=12, height=5)



# ===== Randomized list ===== 

sce.inhib <- sce.backup

# Define the original "Inh_" cell types and their corresponding renaming
inh_celltypes <- c("Inh_5", "Inh_LAMP5_2", "Inh_4", "Inh_LAMP5_1", "Inh_ITC_2", "Inh_ITC_1", "Inh_2", "Inh_SST", "Inh_PVALB", "Inh_1", "Inh_VIP_1", "Inh_VIP_2", "Inh_3")
renaming <- c("PVALB ADAMTS5", "LAMP5 COL14A1", "HTR3A DRD2", "LAMP5 NDNF", "TSHZ1 SEMA3C", "TSHZ1 CALCRL", "VIP NDNF", "SST EPYC", "SST HGF", "VIP NDNF", "VIP ABI3BP", "VIP ABI3BP", "VIP ABI3BP")

# Randomize the renaming
set.seed(1234)  # Setting a seed for reproducibility
random_renaming <- sample(renaming)

# Apply the randomized renaming to the sce.inhib object
for (i in seq_along(inh_celltypes)) {
    colData(sce.inhib)$cell_type[colData(sce.inhib)$cell_type == inh_celltypes[i]] <- random_renaming[i]
}

celltypes_of_interest = sce.inhib[, sce.inhib$cell_type %in% c('TSHZ1 SEMA3C', 'PVALB ADAMTS5', 'VIP ABI3BP',
                                                           'LAMP5 COL14A1', 'HTR3A DRD2', 'LAMP5 NDNF', 
                                                           'TSHZ1 CALCRL', 'VIP NDNF', 'SST EPYC', 
                                                           'SST HGF', 'VIP ABI3BP')]
celltype_matrix = model.matrix(~celltypes_of_interest$cell_type - 1)
colnames(celltype_matrix) = levels(as.factor(celltypes_of_interest$cell_type))

#data(GOhuman)
pdf(here(plot_dir,"Metaneighbor_AUROC_interneurons_shuffled.pdf"))
AUROC_scores_shuffled = MetaNeighbor(dat = celltypes_of_interest ,
                                     experiment_labels = as.numeric(factor(celltypes_of_interest$study_id)),
                                     celltype_labels = celltype_matrix,
                                     genesets = GOhuman,
                                     bplot = TRUE,
                                     fast_version = TRUE)
dev.off()


# plotting using ggplot
long_data <- data.frame(AUROC_scores_shuffled) %>%
    rownames_to_column(var = "GO_term") %>%
    gather(Celltype, AUROC, -GO_term)
head(long_data) 
# GO_term   Celltype AUROC
# 1 GO:0000228 HTR3A.DRD2 0.915
# 2 GO:0000902 HTR3A.DRD2 0.996
# 3 GO:0000988 HTR3A.DRD2 0.964
# 4 GO:0003013 HTR3A.DRD2 0.991
# 5 GO:0003729 HTR3A.DRD2 0.848
# 6 GO:0003735 HTR3A.DRD2 0.657

ggplot(long_data, aes(x = Celltype, y = AUROC)) +
    geom_violin(aes(fill = Celltype), trim = FALSE) +  # Violin plot
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +  # Individual data points
    labs(title = "Violin plot of AUROC scores", y = "AUROC", x = "Celltype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
    coord_cartesian(ylim = c(0, 1))  # Set y-axis limits to 0 and 1


# ============ Directly comparing real vs shuffled data ==============


# Convert AUROC_scores to long format and add a 'Type' column
long_real <- data.frame(AUROC_scores) %>%
    rownames_to_column(var = "GO_term") %>%
    gather(Celltype, AUROC, -GO_term) %>%
    mutate(Type = "Matched")

# Convert AUROC_scores_shuffled to long format and add a 'Type' column
long_shuffled <- data.frame(AUROC_scores_shuffled) %>%
    rownames_to_column(var = "GO_term") %>%
    gather(Celltype, AUROC, -GO_term) %>%
    mutate(Type = "Shuffled")

# Combine the two datasets
combined_data <- bind_rows(long_real, long_shuffled)

# Plot the combined data
ggplot(combined_data, aes(x = Celltype, y = AUROC)) +
    geom_violin(aes(fill = Type), trim = FALSE, position = position_dodge(width = 0.8), alpha=.3) +
    geom_jitter(aes(color = Type), position = position_dodge(width = 0.8), size = 2, alpha = 0.5) +
    stat_summary(aes(group = interaction(Type, Celltype)), fun = mean, geom = "crossbar", color = "inhibck", width=.5, position = position_dodge(width = 0.8)) +
    labs(title = "Comparison of Matched vs Shuffled AUROC scores using GO terms", y = "AUROC", x = "Celltype") +
    scale_fill_manual(values = c("Matched" = "dodgerblue", "Shuffled" = "darkgrey")) +
    scale_color_manual(values = c("Matched" = "dodgerblue", "Shuffled" = "darkgrey")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 1))

ggsave(here(plot_dir, "Metaneighbor_AUROC_interneuron_matched_vs_shuffled.pdf"), width=12, height=5)


# =========== Functional characterization =================

library(GO.db)

# List of GO terms
go_terms <- c("GO:0030674", "GO:0030705", "GO:0043473")

# Retrieve GO term names
go_names <- select(GO.db, keys = as.character(names(GOhuman)), columns = "TERM", keytype = "GOID")

# Print the names
head(go_names)
# GOID                               TERM
# 1 GO:0000228                 nuclear chromosome
# 2 GO:0000902                 cell morphogenesis
# 3 GO:0000988                               <NA>
# 4 GO:0003013         circulatory system process
# 5 GO:0003729                       mRNA binding
# 6 GO:0003735 structural constituent of ribosome


gs_size = sapply(GOhuman, length)
aurocs_df = data.frame(go_term = go_names$TERM, AUROC_scores) 
aurocs_df$average = rowMeans(AUROC_scores)
aurocs_df$n_genes = gs_size[rownames(AUROC_scores)]

head(aurocs_df[order(aurocs_df$average, decreasing = TRUE),],10)
#                                                             go_term HTR3A.DRD2 LAMP5.COL14A1 LAMP5.NDNF
# GO:0000902                                       cell morphogenesis      0.996         0.998      0.995
# GO:0048646 anatomical structure formation involved in morphogenesis      0.993         0.999      0.995
# GO:0005578                                                     <NA>      0.988         0.998      0.993
# GO:0003013                               circulatory system process      0.991         0.995      0.985
# GO:0008283                            cell population proliferation      0.977         0.983      0.975
# GO:0040007                                                   growth      0.978         0.988      0.983
# GO:0008289                                            lipid binding      0.985         0.965      0.971
# GO:0030198                        extracellular matrix organization      0.968         0.996      0.988
# GO:0016757                             glycosyltransferase activity      0.963         0.984      0.978
# GO:0034330                               cell junction organization      0.967         0.974      0.949


small_sets = aurocs_df[aurocs_df$n_genes < 200,] 
head(small_sets[order(small_sets$average, decreasing = TRUE),],10)

#                                                   go_term HTR3A.DRD2 LAMP5.COL14A1 LAMP5.NDNF PVALB.ADAMTS5
# GO:0030674         protein-macromolecule adaptor activity      0.930         0.882      0.924         0.957
# GO:0030705 cytoskeleton-dependent intracellular transport      0.863         0.795      0.916         0.935
# GO:0043473                                   pigmentation      0.824         0.828      0.927         0.897
# GO:0007009                   plasma membrane organization      0.869         0.896      0.829         0.898
# GO:0042393                                histone binding      0.840         0.751      0.872         0.951
# GO:0008565                                           <NA>      0.859         0.840      0.833         0.861
# GO:0007034                             vacuolar transport      0.760         0.714      0.779         0.925
# GO:0032182                 ubiquitin-like protein binding      0.675         0.760      0.757         0.734
# GO:0016779                nucleotidyltransferase activity      0.718         0.756      0.882         0.783
# GO:0005811                                  lipid droplet      0.763         0.748      0.791         0.789