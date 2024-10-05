library("spatialLIBD")
library("patchwork")
library("scran")
library("here")
library("scater")
library("dplyr")
library("tidyr")
library(clusterProfiler)

# Save directories
plot_dir = here("plots", "09_psychiatric_GSEA", "GSEA")
processed_dir <- here("processed-data", "09_psychiatric_GSEA")

# ===== Load human marker genes =====
all_markers <- read.csv(here(processed_dir, "top_100_human_markers_fine_celltype.csv"))

# get top 25 markers genes for each cell type
marker_genes <- all_markers |>
    filter(std.logFC.rank <= 50)


# ====== get Jaffe DEGs =====
jaffe_degs <- read.csv(here(processed_dir,"Jaffe_DEGs.csv"))

# Filter the significant genes (p-value less than 0.05)
PTSD_BLA_genes <- jaffe_degs[jaffe_degs$BLA_PValue_PTSD < 0.05,]
PTSD_BLA_lfc <- jaffe_degs$BLA_logFC_PTSD[jaffe_degs$BLA_PValue_PTSD < 0.05]

MDD_BLA_genes <- jaffe_degs[jaffe_degs$BLA_PValue_MDD < 0.05,]
MDD_BLA_lfc <- jaffe_degs$BLA_logFC_MDD[jaffe_degs$BLA_PValue_MDD < 0.05]


# Prepare the named vector for LRcell
PTSD_degs_vector <- setNames(PTSD_BLA_genes$BLA_PValue_PTSD, PTSD_BLA_genes$Symbol)
length(PTSD_degs_vector)

# ======== Over-representation analysis =========

# === PTSD ===
# Split the genes by cell type
gene_sets <- split(marker_genes$gene, marker_genes$cellType.target)

# Create TERM2GENE data frame
term2gene <- data.frame(
    term = rep(names(gene_sets), lengths(gene_sets)),
    gene = unlist(gene_sets)
)

# Separate DEGs by up-regulated and down-regulated
up_genes <- names(PTSD_degs_vector)[PTSD_BLA_lfc > 0]
length(up_genes)

down_genes <- names(PTSD_degs_vector)[PTSD_BLA_lfc < 0]


# Run enrichment analysis for up-regulated genes
up_enricher_results <- enricher(up_genes, TERM2GENE = term2gene)
# NOTE: no up-enrichment found

# Run enrichment analysis for down-regulated genes
down_enricher_results <- enricher(down_genes, TERM2GENE = term2gene)

# Combine the results into a data frame
down_df <- as.data.frame(down_enricher_results)

# Add direction information
down_df$direction <- "Down-regulated"

png(here(plot_dir, "OverRep_PTSD_downregulated.png"), width = 6, height = 6, units = "in", res = 300)
dotplot(down_enricher_results, showCategory = 20)
dev.off()

# === MDD ===
# Split the genes by cell type
gene_sets <- split(marker_genes$gene, marker_genes$cellType.target)

# Create TERM2GENE data frame
term2gene <- data.frame(
    term = rep(names(gene_sets), lengths(gene_sets)),
    gene = unlist(gene_sets)
)

# Separate DEGs by up-regulated and down-regulated
up_genes <- names(MDD_degs_vector)[MDD_BLA_lfc > 0]
length(up_genes)
# 738

down_genes <- names(MDD_degs_vector)[MDD_BLA_lfc < 0]
length(down_genes)
# 909

# Run enrichment analysis for up-regulated genes
up_enricher_results <- enricher(up_genes, TERM2GENE = term2gene)
# NOTE: no up-enrichment found

# Run enrichment analysis for down-regulated genes
down_enricher_results <- enricher(down_genes, TERM2GENE = term2gene)

# Combine the results into a data frame
up_df <- as.data.frame(up_enricher_results)
down_df <- as.data.frame(down_enricher_results)

# Add direction information
down_df$direction <- "Down-regulated"
up_df$direction <- "Up-regulated"

# Combine the results for visualization
combined_df <- rbind(up_df, down_df)


png(here::here(plot_dir, "OverRep_MDD_combined.png"), width = 5, height = 5, units = "in", res = 300)
ggplot(combined_df, aes(x = direction, y = Description, size = GeneRatio, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  labs(title = "Over-rep with Directionality for MDD DEGs (n=909)", 
       x = "Direction", 
       y = "Cell Type") +
  theme(axis.text.y = element_text(size = 8))  # Adjust text size for better readability
dev.off()



# ====== GSEA =======

PTSD_degs_vector <- sort(PTSD_degs_vector, decreasing = TRUE)  # Assuming PTSD_BLA_lfc is the ranking metric

# Split the genes by cell type (gene sets)
gene_sets <- split(marker_genes$gene, marker_genes$cellType.target)

# Create TERM2GENE data frame
term2gene <- data.frame(
    term = rep(names(gene_sets), lengths(gene_sets)),
    gene = unlist(gene_sets)
)

# Ensure that PTSD_degs_vector is named by gene symbols and ranked by logFC or another metric
# PTSD_degs_vector <- PTSD_degs_vector[order(PTSD_BLA_lfc, decreasing = TRUE)]  # Assuming PTSD_BLA_lfc is the ranking metric

# Run GSEA for the entire ranked gene list (GSEA works with the full ranked list, not just up/down separately)
gsea_results <- GSEA(PTSD_degs_vector, TERM2GENE = term2gene, pvalueCutoff = 0.05)

# Convert GSEA results to a data frame for further manipulation or inspection if needed
gsea_df <- as.data.frame(gsea_results)

# Visualize using dot plot
png(here(plot_dir, "GSEA_PTSD_combined.png"), width = 5, height = 5, units = "in", res = 300)
dotplot(gsea_results, showCategory = 20)
dev.off()
