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
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "08_species_comparisons")


sce.excit <- readRDS(here("processed-data", "07_annotation_and_characterization", "sce_excit_final_subclusters_annotated.rds"))
sce.inhib <- readRDS(here("processed-data", "07_annotation_and_characterization", "sce_inhib_final_subclusters_annotated.rds"))




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

# Define receptor classes and their corresponding genes
serotonin_genes <- c("HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR2A", "HTR2C", "HTR4", "HTR5A", "HTR6", "HTR7")
opioid_genes <- c("OPRM1", "OPRD1", "OPRK1", "OPRL1")
dopamine_genes <- c("DRD1", "DRD2", "DRD3", "DRD4", "DRD5")
norepinephrine_genes <- c("ADRA1A", "ADRA1B", "ADRA2A", "ADRA2C", "ADRB1", "ADRB2")
cannabinoid_genes <- c("CNR1", "CNR2")
#cholinergic_genes <- c("CHRNA3", "CHRNA5", "CHRNA6", "CHRNB3", "CHRNB4")
muscarinic_genes <- c(
  "CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNA6", "CHRNA7", "CHRNA9", "CHRNA10",
  "CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4", 
  "CHRND", "CHRNE", "CHRNG",
  "CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5"
)

# keep only in sce
serotonin_genes <- serotonin_genes[serotonin_genes %in% rownames(sce.excit)]
opioid_genes <- opioid_genes[opioid_genes %in% rownames(sce.excit)]
dopamine_genes <- dopamine_genes[dopamine_genes %in% rownames(sce.excit)]
norepinephrine_genes <- norepinephrine_genes[norepinephrine_genes %in% rownames(sce.excit)]
cannabinoid_genes <- cannabinoid_genes[cannabinoid_genes %in% rownames(sce.excit)]
#cholinergic_genes <- cholinergic_genes[cholinergic_genes %in% rownames(sce.excit)]
muscarinic_genes <- muscarinic_genes[muscarinic_genes %in% rownames(sce.excit)]


# Combine all genes and receptor classes
all_genes <- c(serotonin_genes, opioid_genes, dopamine_genes, norepinephrine_genes, cannabinoid_genes, muscarinic_genes)
receptor_class <- c(
  rep("Serotonin", length(serotonin_genes)),
  rep("Opioid", length(opioid_genes)),
  rep("Dopamine", length(dopamine_genes)),
  rep("Norepinephrine", length(norepinephrine_genes)),
  rep("CBRs", length(cannabinoid_genes)),
  rep("Acetylcholine", length(muscarinic_genes))
)


# Create a data frame to map genes to their receptor class
gene_receptor_map <- data.frame(gene = all_genes, receptor_class = receptor_class)

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

# Merge with receptor class information
percent_data <- merge(percent_expressing_long, gene_receptor_map, by = "gene", all.x = TRUE)

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

# Merge with receptor class information
expr_data_long <- merge(expr_data_long, gene_receptor_map, by = "gene")


# ========== Get subregion info and emerge it all =======
# Extract subcluster information from colData
subcluster_info <- colData(sce.excit_pseudo) %>% as.data.frame()
subcluster_info$fine_celltype <- rownames(subcluster_info)

# Merge the expression data with subcluster information
plot_data <- merge(expr_data_long, subcluster_info, by = "fine_celltype", all.x = TRUE)
plot_data <- merge(plot_data, percent_data, by = c("fine_celltype", "gene"), all.x = TRUE)


# ======== Create the dot plot using ggplot2 ========
png(here(plot_dir, "dotplot_excit_receptor_classes.png"), width = 15, height = 5, units = "in", res = 300)
ggplot(plot_data, aes(x = gene, y = fine_celltype, size = percent_expressing, fill = expression)) +
  geom_point(shape = 21) +  # Dot plot with filled points
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 2.5)) +
  scale_size_continuous(range = c(0, 8)) +  # Adjust the range for the size of the points
  theme_minimal() +
  labs(x = "Gene", y = "Fine Cell Type", fill = "Logcounts", size = "Percent Expressing") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey80", color = "grey50"),
    strip.text = element_text(face = "bold")
  ) +
  facet_grid(subregion_celltype ~ receptor_class.y, scales = "free", space = "free")  # Facet by receptor class and subregion
dev.off()











# ================ Inhibitory ===============


# ======== Faceting around all receptr classes =========

# Define receptor classes and their corresponding genes
serotonin_genes <- c("HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR2A", "HTR2C", "HTR4", "HTR5A", "HTR6", "HTR7")
opioid_genes <- c("OPRM1", "OPRD1", "OPRK1", "OPRL1")
dopamine_genes <- c("DRD1", "DRD2", "DRD3", "DRD4", "DRD5")
norepinephrine_genes <- c("ADRA1A", "ADRA1B", "ADRA2A", "ADRA2C", "ADRB1", "ADRB2")
cannabinoid_genes <- c("CNR1", "CNR2")
#cholinergic_genes <- c("CHRNA3", "CHRNA5", "CHRNA6", "CHRNB3", "CHRNB4")
muscarinic_genes <- c(
  "CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNA6", "CHRNA7", "CHRNA9", "CHRNA10",
  "CHRNB1", "CHRNB2", "CHRNB3", "CHRNB4", 
  "CHRND", "CHRNE", "CHRNG",
  "CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5"
)

# keep only in sce
serotonin_genes <- serotonin_genes[serotonin_genes %in% rownames(sce.inhib)]
opioid_genes <- opioid_genes[opioid_genes %in% rownames(sce.inhib)]
dopamine_genes <- dopamine_genes[dopamine_genes %in% rownames(sce.inhib)]
norepinephrine_genes <- norepinephrine_genes[norepinephrine_genes %in% rownames(sce.inhib)]
cannabinoid_genes <- cannabinoid_genes[cannabinoid_genes %in% rownames(sce.inhib)]
#cholinergic_genes <- cholinergic_genes[cholinergic_genes %in% rownames(sce.inhib)]
muscarinic_genes <- muscarinic_genes[muscarinic_genes %in% rownames(sce.inhib)]


# Combine all genes and receptor classes
all_genes <- c(serotonin_genes, opioid_genes, dopamine_genes, norepinephrine_genes, cannabinoid_genes, muscarinic_genes)
receptor_class <- c(
  rep("Serotonin", length(serotonin_genes)),
  rep("Opioid", length(opioid_genes)),
  rep("Dopamine", length(dopamine_genes)),
  rep("Norepinephrine", length(norepinephrine_genes)),
  rep("CBRs", length(cannabinoid_genes)),
  rep("Acetylcholine", length(muscarinic_genes))
)


# Create a data frame to map genes to their receptor class
gene_receptor_map <- data.frame(gene = all_genes, receptor_class = receptor_class)

# ======== Calculate percentage of cells expressing > 0 counts ========
sce.subset <- sce.inhib[all_genes,]
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

# Merge with receptor class information
percent_data <- merge(percent_expressing_long, gene_receptor_map, by = "gene", all.x = TRUE)

# ========= Pseudobulk for mean logcounts ========
# Aggregate expression data across all defined genes
sce.inhib_pseudo <- aggregateAcrossCells(
  sce.inhib, 
  ids = sce.inhib$fine_celltype, 
  subset_row = all_genes,
  statistics = "mean", 
  use.assay.type = "logcounts"
)

  # Extract expression data from the SingleCellExperiment object
expr_data <- as.data.frame(assay(sce.inhib_pseudo, "logcounts"))

# Add gene and fine cell type information
expr_data$gene <- rownames(expr_data)
expr_data_long <- expr_data %>%
  pivot_longer(cols = -gene, names_to = "fine_celltype", values_to = "expression")

# Merge with receptor class information
expr_data_long <- merge(expr_data_long, gene_receptor_map, by = "gene")


# ========== Get subregion info and emerge it all =======
# Extract subcluster information from colData
subcluster_info <- colData(sce.inhib_pseudo) %>% as.data.frame()
subcluster_info$fine_celltype <- rownames(subcluster_info)

# Merge the expression data with subcluster information
plot_data <- merge(expr_data_long, subcluster_info, by = "fine_celltype", all.x = TRUE)
plot_data <- merge(plot_data, percent_data, by = c("fine_celltype", "gene"), all.x = TRUE)



# Create a dot plot using ggplot2
png(here(plot_dir, "dotplot_inhib_receptor_classes.png"), width = 15, height = 5, units = "in", res = 300)
ggplot(plot_data, aes(x = gene, y = fine_celltype, size=percent_expressing, fill = expression)) +
  geom_point(shape = 21) +  # Dot plot with filled points
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 2.5), , oob = scales::squish) +
  scale_size_continuous(range = c(0, 7)) +  # Adjust the range for the size of the points
  theme_minimal() +
  labs(x = "Gene", y = "Fine Cell Type", fill = "Logcounts", size = "Percent Expressing") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey80", color = "grey50"),
    strip.text = element_text(face = "bold")
  ) +
  facet_grid( ~receptor_class.y, scales = "free", space = "free")  # Facet by receptor class and subregion
dev.off()



