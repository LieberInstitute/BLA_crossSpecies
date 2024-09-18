library(dplyr)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(scater)
library(here)

# directories
processed_dir = here("processed-data/07_annotation_and_characterization")
plot_dir = here("plots", "08_species_comparisons")


sce.excit <- readRDS(here(processed_dir, "sce_excit_final_subclusters_annotated.rds"))
sce.inhib <- readRDS(here(processed_dir, "sce_inhib_final_subclusters_annotated.rds"))



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


# markers.baboon$t.stat <- markers.baboon$std.logFC * sqrt(ncol(sce.excit.baboon))
# markers.human$t.stat <- markers.human$std.logFC * sqrt(ncol(sce.excit.human))
# markers.macaque$t.stat <- markers.macaque$std.logFC * sqrt(ncol(sce.excit.macaque))



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
                            std.logFC = markers.human_PEX5L_MYRIP$std.logFC)

markers.macaque.df <- data.frame(gene = markers.macaque_PEX5L_MYRIP$gene,
                              std.logFC = markers.macaque_PEX5L_MYRIP$std.logFC)

markers.baboon.df <- data.frame(gene = markers.baboon_PEX5L_MYRIP$gene,
                             std.logFC = markers.baboon_PEX5L_MYRIP$std.logFC)

# Rename the columns for clarity before merging
colnames(markers.human.df)[2] <- "human_std.logFC"
colnames(markers.macaque.df)[2] <- "macaque_std.logFC"
colnames(markers.baboon.df)[2] <- "baboon_std.logFC"

# Merge the dataframes on the gene column
merged.df <- merge(markers.human.df, markers.macaque.df, by = "gene", all = TRUE)
merged.df <- merge(merged.df, markers.baboon.df, by = "gene", all = TRUE)

# Display the merged dataframe
print(merged.df)

# get correlation
cor_val <- cor(merged.df$macaque_std.logFC, merged.df$human_std.logFC)

# plot using ggplot

# set seed
set.seed(1234)
png(here(plot_dir, "correlation_std.logFC_LA.png"), width = 3, height = 3, units = "in", res = 300)
ggplot(merged.df, aes(x = human_std.logFC, y = macaque_std.logFC )) +
    geom_point(size = 0.5) +  # Default points in black
    theme_minimal() +
    labs(x = "Human LA std.logFC", y = "Macaque LA std.logFC") +
    ggtitle(paste0("r = ", print(plyr::round_any(cor_val, .001)))) +
    xlim(-3, 3) +
    ylim(-3, 3) +
    # Add horizontal and vertical lines
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme(panel.border = element_rect(colour = "black", fill = NA)) + 
    # Add red points for GULP1 and SATB1
    geom_point(data = merged.df[merged.df$gene %in% c("GULP1", "SATB1"),], 
               aes(x = human_std.logFC, y = macaque_std.logFC), 
               color = "#00BFC4", size = 1.5) +  # Make specific points red and slightly larger
    # Add red labels for GULP1 and SATB1
    geom_text_repel(data = merged.df[merged.df$gene %in% c("GULP1", "SATB1"),], 
                    aes(label = gene), 
                    color = "#00BFC4", 
                    box.padding = .2,
                    fontface = "bold")
dev.off()



# ======= Basal  =========
celltype.of.interest <- "BA"

markers.human_PEX5L_MYRIP <- markers.human %>%
    filter(cellType.target == celltype.of.interest)

markers.macaque_PEX5L_MYRIP <- markers.macaque %>%
    filter(cellType.target == celltype.of.interest)

markers.baboon_PEX5L_MYRIP <- markers.baboon %>%
    filter(cellType.target == celltype.of.interest)

# make dataframe of gene, human_std.logFC
markers.human.df <- data.frame(gene = markers.human_PEX5L_MYRIP$gene,
                               std.logFC = markers.human_PEX5L_MYRIP$std.logFC)

markers.macaque.df <- data.frame(gene = markers.macaque_PEX5L_MYRIP$gene,
                                 std.logFC = markers.macaque_PEX5L_MYRIP$std.logFC)

markers.baboon.df <- data.frame(gene = markers.baboon_PEX5L_MYRIP$gene,
                                std.logFC = markers.baboon_PEX5L_MYRIP$std.logFC)

# Rename the columns for clarity before merging
colnames(markers.human.df)[2] <- "human_std.logFC"
colnames(markers.macaque.df)[2] <- "macaque_std.logFC"
colnames(markers.baboon.df)[2] <- "baboon_std.logFC"

# Merge the dataframes on the gene column
merged.df <- merge(markers.human.df, markers.macaque.df, by = "gene", all = TRUE)
merged.df <- merge(merged.df, markers.baboon.df, by = "gene", all = TRUE)

# Display the merged dataframe
print(merged.df)

# get correlation
cor_val <- cor(merged.df$macaque_std.logFC, merged.df$human_std.logFC)

# plot using ggplot

png(here(plot_dir, "correlation_std.logFC_BA.png"), width = 3, height = 3, units = "in", res = 300)
ggplot(merged.df, aes(x = human_std.logFC, y = macaque_std.logFC )) +
    geom_point(size=.5) +
    theme_minimal() +
    labs(x = "Human BA std.logFC", y = "Macaque BA std.logFC") +
    ggtitle(paste0("r = ", print(plyr::round_any(cor_val, .001)))) +
    xlim(-3, 3) +
    ylim(-3, 3) +
    geom_text_repel(data = merged.df[order(merged.df$macaque_std.logFC, decreasing = TRUE),][1:10,],
                    aes(label = gene), box.padding = 0.2) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme(panel.border = element_rect(colour = "black", fill=NA)) 
dev.off()


# set seed
set.seed(18995)
png(here(plot_dir, "correlation_std.logFC_BA.png"), width = 3, height = 3, units = "in", res = 300)
ggplot(merged.df, aes(x = human_std.logFC, y = macaque_std.logFC )) +
    geom_point(size = 0.5) +  # Default points in black
    theme_minimal() +
    labs(x = "Human BA std.logFC", y = "Macaque BA std.logFC") +
    ggtitle(paste0("r = ", print(plyr::round_any(cor_val, .001)))) +
    xlim(-3, 3) +
    ylim(-3, 3) +
    # Add horizontal and vertical lines
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme(panel.border = element_rect(colour = "black", fill = NA)) + 
    # Add red points for GULP1 and SATB1
    geom_point(data = merged.df[merged.df$gene %in% c("COL25A1",  "PEX5L"),], 
               aes(x = human_std.logFC, y = macaque_std.logFC), 
               color = "#F8766D", size = 1.5) +  # Make specific points red and slightly larger
    # Add red labels for GULP1 and SATB1
    geom_text_repel(data = merged.df[merged.df$gene %in% c("COL25A1",  "PEX5L"),], 
                    aes(label = gene), 
                    box.padding = 0.2,
                    color="#F8766D",
                    fontface = "bold"
                    )
dev.off()




# ======= Accessort Basal  =========
celltype.of.interest <- "aBA"

markers.human_PEX5L_MYRIP <- markers.human %>%
    filter(cellType.target == celltype.of.interest)

markers.macaque_PEX5L_MYRIP <- markers.macaque %>%
    filter(cellType.target == celltype.of.interest)

markers.baboon_PEX5L_MYRIP <- markers.baboon %>%
    filter(cellType.target == celltype.of.interest)

# make dataframe of gene, human_std.logFC
markers.human.df <- data.frame(gene = markers.human_PEX5L_MYRIP$gene,
                               std.logFC = markers.human_PEX5L_MYRIP$std.logFC)

markers.macaque.df <- data.frame(gene = markers.macaque_PEX5L_MYRIP$gene,
                                 std.logFC = markers.macaque_PEX5L_MYRIP$std.logFC)

markers.baboon.df <- data.frame(gene = markers.baboon_PEX5L_MYRIP$gene,
                                std.logFC = markers.baboon_PEX5L_MYRIP$std.logFC)

# Rename the columns for clarity before merging
colnames(markers.human.df)[2] <- "human_std.logFC"
colnames(markers.macaque.df)[2] <- "macaque_std.logFC"
colnames(markers.baboon.df)[2] <- "baboon_std.logFC"

# Merge the dataframes on the gene column
merged.df <- merge(markers.human.df, markers.macaque.df, by = "gene", all = TRUE)
merged.df <- merge(merged.df, markers.baboon.df, by = "gene", all = TRUE)

# Display the merged dataframe
print(merged.df)

# get correlation
cor_val <- cor(merged.df$macaque_std.logFC, merged.df$human_std.logFC)

# plot using ggplot

# png(here(plot_dir, "correlation_std.logFC_aBA.png"), width = 3, height = 3, units = "in", res = 300)
# ggplot(merged.df, aes(x = human_std.logFC, y = macaque_std.logFC )) +
#     geom_point(size=.5) +
#     theme_minimal() +
#     labs(x = "Human aBA std.logFCistic", y = "Macaque aBA std.logFCistic") +
#     ggtitle(paste0("r = ", print(plyr::round_any(cor_val, .001)))) +
#     xlim(-3, 3) +
#     ylim(-3, 3) +
#     geom_text_repel(data = merged.df[order(merged.df$macaque_std.logFC, decreasing = TRUE),][1:10,],
#                   aes(label = gene), box.padding = 0.5) +
#     geom_hline(yintercept = 0, linetype = "dotted") +
#     geom_vline(xintercept = 0, linetype = "dotted") +
#     theme(panel.border = element_rect(colour = "black", fill=NA)) 
# dev.off()



set.seed(18995)
png(here(plot_dir, "correlation_std.logFC_aBA.png"), width = 3, height = 3, units = "in", res = 300)
ggplot(merged.df, aes(x = human_std.logFC, y = macaque_std.logFC )) +
    geom_point(size = 0.5) +  # Default points in black
    theme_minimal() +
    labs(x = "Human aBA std.logFC", y = "Macaque aBA std.logFC") +
    ggtitle(paste0("r = ", print(plyr::round_any(cor_val, .001)))) +
    xlim(-3, 3) +
    ylim(-3, 3) +
    # Add horizontal and vertical lines
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme(panel.border = element_rect(colour = "black", fill = NA)) + 
    # Add red points for GULP1 and SATB1
    geom_point(data = merged.df[merged.df$gene %in% c("ESR1"),], 
               aes(x = human_std.logFC, y = macaque_std.logFC), 
               color = "#C77CFF", size = 1.5) +  # Make specific points red and slightly larger
    # Add red labels for GULP1 and SATB1
    geom_text_repel(data = merged.df[merged.df$gene %in% c("ESR1"),], 
                    aes(label = gene), 
                    box.padding = 0.2,
                    color="#C77CFF",
                    fontface = "bold"
                    )
dev.off()
