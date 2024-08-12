library(dplyr)
library(here)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(scater)

# directories
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "08_species_comparisons")


sce.excit <- readRDS(here("processed-data", "sce_excit_final_subclusters_annotated.rds"))
sce.inhib <- readRDS(here("processed-data", "sce_inhib_final_subclusters_annotated.rds"))



# ========= Markers genes per species =========
library(DeconvoBuddies)

# combine putative LA, BA, and aBA clusters
LA_clusters <- c("GULP1_TRHDE", "ZBTB20_SLC4A4")
BA_clusters <- c("PEX5L_MYRIP", "MEIS2_COL25A1")
aBA_clusters <- c("ESR1_ADRA1A", "GRIK3_TNS3")

sce.excit$subregion_celltype <- sce.excit$fine_celltype
sce.excit$subregion_celltype[sce.excit$fine_celltype %in% LA_clusters] <- "LA"
sce.excit$subregion_celltype[sce.excit$fine_celltype %in% BA_clusters] <- "BA"
sce.excit$subregion_celltype[sce.excit$fine_celltype %in% aBA_clusters] <- "aBA"

#subset species
sce.excit.human <- sce.excit[,which(colData(sce.excit)$species == "human")]
sce.excit.macaque <- sce.excit[,which(colData(sce.excit)$species == "macaque")]
sce.excit.baboon <- sce.excit[,which(colData(sce.excit)$species == "baboon")]


# find markers
markers.human <- findMarkers_1vAll(sce.excit.human,
                                   assay_name = "logcounts",
                                   cellType_col = "subregion_celltype",
                                   add_symbol = FALSE,
                                   mod = "~Sample",
                                   verbose = TRUE)

markers.macaque <- findMarkers_1vAll(sce.excit.macaque,
                                     assay_name = "logcounts",
                                     cellType_col = "subregion_celltype",
                                     add_symbol = FALSE,
                                     mod = "~Subject",
                                     verbose = TRUE)

markers.baboon <- findMarkers_1vAll(sce.excit.baboon,
                                    assay_name = "logcounts",
                                    cellType_col = "subregion_celltype",
                                    add_symbol = FALSE,
                                    mod = "~Sample",
                                    verbose = TRUE)


# MORE NOTES:
#
# I forgot to define putative BA vs LA vs aBA celltypes
# This would be a more interesting comparison of marker genes. 




# convert to t-statistic

round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}


markers.baboon$t.stat <- markers.baboon$std.logFC * sqrt(ncol(sce.excit.baboon))
markers.human$t.stat <- markers.human$std.logFC * sqrt(ncol(sce.excit.human))
markers.macaque$t.stat <- markers.macaque$std.logFC * sqrt(ncol(sce.excit.macaque))



# ======= Lateral =========
celltype.of.interest <- "LA"

markers.human_PEX5L_MYRIP <- markers.human %>%
    filter(cellType.target == celltype.of.interest)

markers.macaque_PEX5L_MYRIP <- markers.macaque %>%
    filter(cellType.target == celltype.of.interest)

markers.baboon_PEX5L_MYRIP <- markers.baboon %>%
    filter(cellType.target == celltype.of.interest)

# make dataframe of gene, human_t.stat
markers.human.df <- data.frame(gene = markers.human_PEX5L_MYRIP$gene,
                            t.stat = markers.human_PEX5L_MYRIP$t.stat)

markers.macaque.df <- data.frame(gene = markers.macaque_PEX5L_MYRIP$gene,
                              t.stat = markers.macaque_PEX5L_MYRIP$t.stat)

markers.baboon.df <- data.frame(gene = markers.baboon_PEX5L_MYRIP$gene,
                             t.stat = markers.baboon_PEX5L_MYRIP$t.stat)

# Rename the columns for clarity before merging
colnames(markers.human.df)[2] <- "human_t.stat"
colnames(markers.macaque.df)[2] <- "macaque_t.stat"
colnames(markers.baboon.df)[2] <- "baboon_t.stat"

# Merge the dataframes on the gene column
merged.df <- merge(markers.human.df, markers.macaque.df, by = "gene", all = TRUE)
merged.df <- merge(merged.df, markers.baboon.df, by = "gene", all = TRUE)

# Display the merged dataframe
print(merged.df)

# get correlation
cor_val <- cor(merged.df$macaque_t.stat, merged.df$human_t.stat)

# plot using ggplot

png(here(plot_dir, "correlation_t.stat_LA.png"), width = 3, height = 3, units = "in", res = 300)
ggplot(merged.df, aes(x = human_t.stat, y = macaque_t.stat )) +
    geom_point(size=.5) +
    theme_minimal() +
    labs(x = "Human LA t-statistic", y = "Macaque LA t-statistic") +
    ggtitle(paste0("r = ", print(round_any(cor_val, .001)))) +
    xlim(-400, 400) +
    ylim(-400, 400) +
    #geom_text_repel(data = merged.df[order(merged.df$macaque_t.stat, decreasing = TRUE),][1:5,],
                    #aes(label = gene), box.padding = 0.2) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme(panel.border = element_rect(colour = "black", fill=NA))+ 
    # show GULP1 and SATB1
    geom_text_repel(data = merged.df[merged.df$gene %in% c("GULP1", "SATB1"),], aes(label = gene), box.padding = 0.2)
dev.off()




# ======= Basal  =========
celltype.of.interest <- "BA"

markers.human_PEX5L_MYRIP <- markers.human %>%
    filter(cellType.target == celltype.of.interest)

markers.macaque_PEX5L_MYRIP <- markers.macaque %>%
    filter(cellType.target == celltype.of.interest)

markers.baboon_PEX5L_MYRIP <- markers.baboon %>%
    filter(cellType.target == celltype.of.interest)

# make dataframe of gene, human_t.stat
markers.human.df <- data.frame(gene = markers.human_PEX5L_MYRIP$gene,
                               t.stat = markers.human_PEX5L_MYRIP$t.stat)

markers.macaque.df <- data.frame(gene = markers.macaque_PEX5L_MYRIP$gene,
                                 t.stat = markers.macaque_PEX5L_MYRIP$t.stat)

markers.baboon.df <- data.frame(gene = markers.baboon_PEX5L_MYRIP$gene,
                                t.stat = markers.baboon_PEX5L_MYRIP$t.stat)

# Rename the columns for clarity before merging
colnames(markers.human.df)[2] <- "human_t.stat"
colnames(markers.macaque.df)[2] <- "macaque_t.stat"
colnames(markers.baboon.df)[2] <- "baboon_t.stat"

# Merge the dataframes on the gene column
merged.df <- merge(markers.human.df, markers.macaque.df, by = "gene", all = TRUE)
merged.df <- merge(merged.df, markers.baboon.df, by = "gene", all = TRUE)

# Display the merged dataframe
print(merged.df)

# get correlation
cor_val <- cor(merged.df$macaque_t.stat, merged.df$human_t.stat)

# plot using ggplot

png(here(plot_dir, "correlation_t.stat_BA.png"), width = 3, height = 3, units = "in", res = 300)
ggplot(merged.df, aes(x = human_t.stat, y = macaque_t.stat )) +
    geom_point(size=.5) +
    theme_minimal() +
    labs(x = "Human BA t-statistic", y = "Macaque BA t-statistic") +
    ggtitle(paste0("r = ", print(round_any(cor_val, .001)))) +
    xlim(-400, 400) +
    ylim(-400, 400) +
    geom_text_repel(data = merged.df[order(merged.df$macaque_t.stat, decreasing = TRUE),][1:10,],
                    aes(label = gene), box.padding = 0.2) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme(panel.border = element_rect(colour = "black", fill=NA)) 
dev.off()



# ======= Accessort Basal  =========
celltype.of.interest <- "aBA"

markers.human_PEX5L_MYRIP <- markers.human %>%
    filter(cellType.target == celltype.of.interest)

markers.macaque_PEX5L_MYRIP <- markers.macaque %>%
    filter(cellType.target == celltype.of.interest)

markers.baboon_PEX5L_MYRIP <- markers.baboon %>%
    filter(cellType.target == celltype.of.interest)

# make dataframe of gene, human_t.stat
markers.human.df <- data.frame(gene = markers.human_PEX5L_MYRIP$gene,
                               t.stat = markers.human_PEX5L_MYRIP$t.stat)

markers.macaque.df <- data.frame(gene = markers.macaque_PEX5L_MYRIP$gene,
                                 t.stat = markers.macaque_PEX5L_MYRIP$t.stat)

markers.baboon.df <- data.frame(gene = markers.baboon_PEX5L_MYRIP$gene,
                                t.stat = markers.baboon_PEX5L_MYRIP$t.stat)

# Rename the columns for clarity before merging
colnames(markers.human.df)[2] <- "human_t.stat"
colnames(markers.macaque.df)[2] <- "macaque_t.stat"
colnames(markers.baboon.df)[2] <- "baboon_t.stat"

# Merge the dataframes on the gene column
merged.df <- merge(markers.human.df, markers.macaque.df, by = "gene", all = TRUE)
merged.df <- merge(merged.df, markers.baboon.df, by = "gene", all = TRUE)

# Display the merged dataframe
print(merged.df)

# get correlation
cor_val <- cor(merged.df$macaque_t.stat, merged.df$human_t.stat)

# plot using ggplot

png(here(plot_dir, "correlation_t.stat_aBA.png"), width = 3, height = 3, units = "in", res = 300)
ggplot(merged.df, aes(x = human_t.stat, y = macaque_t.stat )) +
    geom_point(size=.5) +
    theme_minimal() +
    labs(x = "Human aBA t-statistic", y = "Macaque aBA t-statistic") +
    ggtitle(paste0("r = ", print(round_any(cor_val, .001)))) +
    xlim(-400, 400) +
    ylim(-400, 400) +
    geom_text_repel(data = merged.df[order(merged.df$macaque_t.stat, decreasing = TRUE),][1:10,],
                  aes(label = gene), box.padding = 0.2) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme(panel.border = element_rect(colour = "black", fill=NA)) 
dev.off()




genes <- c(
    "GRM1", "GRM5", "GRM2", "GRM3", "GRM4", "GRM6", "GRM7", "GRM8",
    "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5",
    "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B",
    "HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR2A", "HTR2C", "HTR4", "HTR5A", "HTR6", "HTR7",
    "DRD1", "DRD2", "DRD3", "DRD4", "DRD5",
    "CHRNA3", "CHRNA5", "CHRNA6", "CHRNB3", "CHRNB4",
    "ADRA1A", "ADRA1B", "ADRA2A", "ADRA2C", "ADRB1", "ADRB2", "AVPR1A", "AVPR1B", "AVPR2", "CNR1", "CNR2", 
    "DRD2", "DRD3", "F2RL1", "F2RL2", "F2RL3", "GHRHR", "HCRTR2", "KISS1R", "MC2R", "MC3R", "MC4R", "MC5R", 
    "NPBWR1", "NPBWR2", "NPY1R", "NPY2R", "NPY4R", "NPY5R", "OPRD1", "OPRK1", "OPRM1", "OXTR", "PRLHR", 
    "PTGER1", "PTGER2", "PTGER3", "PTGER4", "PTH1R", "PTH2R", "SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5", 
    "TACR1", "TACR2", "TACR3", "TRHR", "V1AR", "V1BR"
)

# drop any not in sce.excit
genes <- genes[genes %in% rownames(sce.excit)]

# drop any duplicate
genes <- unique(genes)


png(here(plot_dir, "dotplot_excitatory_receptrs_and_such.png"), width = 20, height = 5, units = "in", res = 300)
plotDots(sce.excit, features=genes, group="fine_celltype") +
    coord_flip() +
    scale_color_gradient(low = "white", high = "red") +
    # x label horizontal
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()