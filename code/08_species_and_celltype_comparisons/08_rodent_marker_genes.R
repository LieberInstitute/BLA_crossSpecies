library(here)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(dplyr)
library(scuttle)
library(tidyr)


# directories
processed_dir = here("processed-data", "08_species_comparisons", "06_rodent_markers")
plot_dir = here("plots", "08_species_comparisons", "06_rodent_markers")


sce.excit <- readRDS(here("processed-data", "07_annotation_and_characterization", "sce_excit_final_subclusters_annotated.rds"))
# read csv for mose markers

mouse_markers <- read.csv(here(processed_dir, "mousegenes.txt"), header=FALSE)
mouse_markers

#         V1
# 1    GULP1
# 2  COL25A1
# 3    CPLX1
# 4    RSPO2
# 5  CAMK2N1
# 6  COL12A1
# 7    CPNE8
# 8    NEGR1
# 9     RORB
# 10     GRP
# 11   MEF2C
# 12   SYNPR
# 13   TSHZ2
# 14 ADAMTS2
# 15    ETV1
# 16  GRIN3A
# 17   OPRK1
# 18 SLC24A3
# 19  SEMA3E
# 20   SLIT2
# 21    SOX5


png(here(plot_dir, "mouse_markers.png"), width = 6, height = 12, units = "in", res = 300)
plotExpression(sce.excit,
              features=mouse_markers$V1,
              x="fine_celltype",
              colour_by="fine_celltype",
              ncol=2
              )
dev.off()


png(here(plot_dir, "mouse_markers_dots.png"), width = 10, height = 7, units = "in", res = 300)
plotDots(sce.excit,
         features=mouse_markers$V1,
         group="fine_celltype",
         color=c("white", "red"),
         scale=TRUE
         ) +
         coord_flip() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# ======= Making a better dotplot =======



# ================ Excitatory ===============

# combine putative LA, BA, and aBA clusters
LA_clusters <- c("GULP1_TRHDE", "ZBTB20_SLC4A4")
BA_clusters <- c("PEX5L_MYRIP", "MEIS2_COL25A1")
aBA_clusters <- c("ESR1_ADRA1A", "GRIK3_TNS3")

sce.excit$subregion_celltype <- "Other"
sce.excit$subregion_celltype[sce.excit$fine_celltype %in% LA_clusters] <- "LA"
sce.excit$subregion_celltype[sce.excit$fine_celltype %in% BA_clusters] <- "BA"
sce.excit$subregion_celltype[sce.excit$fine_celltype %in% aBA_clusters] <- "aBA"
unique(sce.excit$subregion_celltype)
# [1] "aBA"   "BA"    "Other" "LA"   


# ======== Faceting around all receptr classes =========



# Define marker classes and their corresponding genes
primate_markers <- c("GULP1", "SATB1", "COL25A1", "PEX5L", "ESR1")
# drop GULP1 and COL25A1 from mouse_markers
rodent_markers <- mouse_markers[-c(1, 2),]
rodent_markers <- c(rodent_markers, "PPP1R1B")

# Combine all genes and marker classes
all_genes <- c(primate_markers, rodent_markers)
marker_class <- c(
  rep("Primate Markers", length(primate_markers)),
  rep("Rodent Markers", length(rodent_markers))
)


# Create a data frame to map genes to their marker class
gene_marker_map <- data.frame(gene = all_genes, marker_class = marker_class)

# ======== Calculate percentage of cells expressing > 0 counts ========
sce.subset <- sce.excit[all_genes,]
percent_expressing <- apply(assay(sce.subset, "counts") > 0, 1, function(x) {
  tapply(x, sce.subset$fine_celltype, mean) * 100
})

# Convert to a long format data frame
percent_expressing_df <- as.data.frame(percent_expressing)
percent_expressing_df$gene <- rownames(percent_expressing_df)
percent_expressing_long <- percent_expressing_df %>%
  pivot_longer(cols = -gene, names_to = "fine_celltype", values_to = "percent_expressing")

# change gene column name to "fine_celltype", and visa vers
percent_expressing_long <- percent_expressing_long %>% rename(fine_celltype = gene, gene = fine_celltype)

# Merge with marker class information
percent_data <- merge(percent_expressing_long, gene_marker_map, by = "gene", all.x = TRUE)

# ========= Pseudobulk for mean logcounts ========
# Aggregate expression data across all defined genes
sce.excit_pseudo <- aggregateAcrossCells(
  sce.excit, 
  ids = sce.excit$fine_celltype, 
  subset_row = all_genes,
  statistics = "mean", 
  use.assay.type = "logcounts"
)

  # Extract expression data from the SingleCellExperiment object
expr_data <- as.data.frame(assay(sce.excit_pseudo, "logcounts"))

# Add gene and fine cell type information
expr_data$gene <- rownames(expr_data)
expr_data_long <- expr_data %>%
  pivot_longer(cols = -gene, names_to = "fine_celltype", values_to = "expression")

# Merge with marker class information
expr_data_long <- merge(expr_data_long, gene_marker_map, by = "gene")


# ========== Get subregion info and emerge it all =======
# Extract subcluster information from colData
subcluster_info <- colData(sce.excit_pseudo) %>% as.data.frame()
subcluster_info$fine_celltype <- rownames(subcluster_info)

# Merge the expression data with subcluster information
plot_data <- merge(expr_data_long, subcluster_info, by = "fine_celltype", all.x = TRUE)
plot_data <- merge(plot_data, percent_data, by = c("fine_celltype", "gene"), all.x = TRUE)


# ======== Create the dot plot using ggplot2 ========
png(here(plot_dir, "dotplot_excit_marker_classes.png"), width = 15, height = 5, units = "in", res = 300)
ggplot(plot_data, aes(x = gene, y = fine_celltype, size = percent_expressing, fill = expression)) +
  geom_point(shape = 21) +  # Dot plot with filled points
  scale_fill_gradient(low = "white", high = "red") +
  scale_size_continuous(range = c(0, 8)) +  # Adjust the range for the size of the points
  theme_minimal() +
  labs(x = "Gene", y = "Fine Cell Type", fill = "Logcounts", size = "Percent Expressing") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey80", color = "grey50"),
    strip.text = element_text(face = "bold", size = 14)
  ) +
  facet_grid(subregion_celltype ~ marker_class.y, scales = "free", space = "free")  # Facet by marker class and subregion
dev.off()
