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

# directories
raw_dir = here("raw-data")
samples_dir = here(raw_dir, "samples")
sampleinfo_dir = here(raw_dir, "sampleinfo")
processed_dir = here("processed-data")
plot_dir = here("plots", "00_costa_rds")

Amy.700 <- readRDS(here(processed_dir, "00_costa_rds", "amygdalaCNScaled700.rds"))
Amy.700
sce.costa <- as.SingleCellExperiment(Amy.700)
# An object of class Seurat 
# 21353 features across 120794 samples within 2 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

Bab.200 <-readRDS(here(processed_dir, "00_costa_rds", "AllBaboon_Amygdala200_04.rds"))
Bab.200
sce.bab <- as.SingleCellExperiment(Bab.200)

human.yu <- readRDS(here(processed_dir, "00_costa_rds", "GSE195445_Human_obj.rda"))
human.yu
sce.yu <- as.SingleCellExperiment(human.yu)

# remove "Human" form the beginning of $ident
levels(sce.yu$ident) <- gsub("Human_", "", levels(sce.yu$ident))

load(here(processed_dir, "00_costa_rds", "sce_annotated.rda"))
sce.totty <- sce

# ======= Pearson correlation within Macaque celltypes using all genes =======

library(reshape2)

# 1. Extract expression data
Mac_exprs_data <- logcounts(sce.costa)  # or use `assay(sce, "normalized")` or similar depending on your data
Bab_exprs_data <- logcounts(sce.bab)  # or use `assay(sce, "normalized")` or similar depending on your data
Totty_exprs_data <- logcounts(sce.totty)  # or use `assay(sce, "normalized")` or similar depending on your data
Yu_exprs_data <- logcounts(sce.yu)  # or use `assay(sce, "normalized")` or similar depending on your data

# 2. Identify common genes
common_genes <- Reduce(intersect, list(rownames(Mac_exprs_data), 
                                       rownames(Totty_exprs_data), 
                                       rownames(Yu_exprs_data),
                                       rownames(Bab_exprs_data)))

# 3. Subset each dataset
Mac_exprs_data_common <- Mac_exprs_data[common_genes, ]
Totty_exprs_data_common <- Totty_exprs_data[common_genes, ]
Yu_exprs_data_common <- Yu_exprs_data[common_genes, ]
Bab_exprs_data_common <- Bab_exprs_data[common_genes, ]

# scale if you want to
Mac_exprs_data_common <- scale(Mac_exprs_data_common)
Totty_exprs_data_common <- scale(Totty_exprs_data_common)
Yu_exprs_data_common <- scale(Yu_exprs_data_common )
Bab_exprs_data_common <- scale(Bab_exprs_data_common )


# 4. Group by clusters and compute mean expression
Mac_clusters <- unique(colData(sce.costa)$seurat_clusters)
Totty_clusters <- unique(colData(sce.totty)$annotation)
Yu_clusters <- unique(colData(sce.yu)$ident)
Bab_clusters <- unique(colData(sce.bab)$ident)



#  function to get mean expression for each cluster
compute_mean_exprs <- function(exprs_data, sce, cluster_col, cluster) {
  subset_data <- exprs_data[, colData(sce)[[cluster_col]] == cluster]
  return(rowMeans(subset_data))
}

# Now use the function for each dataset:

# For Mac data
Mac_mean_exprs_list <- lapply(Mac_clusters, function(cluster) {
  compute_mean_exprs(Mac_exprs_data_common, sce.costa, "seurat_clusters", cluster)
})

# For Totty data
Totty_mean_exprs_list <- lapply(Totty_clusters, function(cluster) {
  compute_mean_exprs(Totty_exprs_data_common, sce.totty, "annotation", cluster)
})

# For Yu data
Yu_mean_exprs_list <- lapply(Yu_clusters, function(cluster) {
  compute_mean_exprs(Yu_exprs_data_common, sce.yu, "ident", cluster)
})

# For Baboon data
Bab_mean_exprs_list <- lapply(Bab_clusters, function(cluster) {
  compute_mean_exprs(Bab_exprs_data_common, sce.bab, "ident", cluster)
})


# Convert the list of mean expressions to a matrix
Mac_mean_exprs_matrix <- do.call(cbind, Mac_mean_exprs_list)
colnames(Mac_mean_exprs_matrix) <- Mac_clusters

Totty_mean_exprs_matrix <- do.call(cbind, Totty_mean_exprs_list)
colnames(Totty_mean_exprs_matrix) <- Totty_clusters

Yu_mean_exprs_matrix <- do.call(cbind, Yu_mean_exprs_list)
colnames(Yu_mean_exprs_matrix) <- Yu_clusters

Bab_mean_exprs_matrix <- do.call(cbind, Bab_mean_exprs_list)
colnames(Bab_mean_exprs_matrix) <- Bab_clusters




# 5. Calculate Pearson correlation
Mac_Totty_cor_matrix <- cor(Mac_mean_exprs_matrix, Totty_mean_exprs_matrix)
Mac_Yu_cor_matrix <- cor(Mac_mean_exprs_matrix, Yu_mean_exprs_matrix)
Totty_Yu_cor_matrix <- cor(Totty_mean_exprs_matrix, Yu_mean_exprs_matrix)

Bab_Mac_cor_matrix <- cor(Bab_mean_exprs_matrix, Mac_mean_exprs_matrix)
Bab_Totty_cor_matrix <- cor(Bab_mean_exprs_matrix, Totty_mean_exprs_matrix)
Bab_Yu_cor_matrix <- cor(Bab_mean_exprs_matrix, Yu_mean_exprs_matrix)


# 5. Visualize correlation


# ========== Get annotations ==========

# === Totty celltypes ===
sce.totty$celltype <- factor(sce.totty$celltype)
levels(sce.totty$celltype) <- gsub("-", "_", levels(sce.totty$celltype))

# Create a dataframe from the SCE object
cluster_annotations_df <- data.frame(
  annotation = sce.totty$annotation,
  Totty_celltype = sce.totty$celltype
)

# Remove duplicate rows to get unique mappings
Totty_cluster_annotations_df <- unique(cluster_annotations_df)

# Set the 'annotation' column as rownames
rownames(Totty_cluster_annotations_df) <- Totty_cluster_annotations_df$annotation
Totty_cluster_annotations_df$annotation <- NULL

# View the dataframe
head(Totty_cluster_annotations_df)
# Totty_celltype
# Inh_LAMP5_2          Inhib
# Exc_09               Excit
# Inh_3                Inhib
# Exc_03               Excit
# Inh_LAMP5_1          Inhib
# Inh_ITC_1            Inhib



# === Yu celltypes ===

# Create a dataframe from the SCE object
cluster_annotations_df <- data.frame(
  annotation = sce.yu$ident,
  Yu_celltype = sce.yu$celltype
)

# Remove duplicate rows to get unique mappings
Yu_cluster_annotations_df <- unique(cluster_annotations_df)

# Set the 'annotation' column as rownames
rownames(Yu_cluster_annotations_df) <- Yu_cluster_annotations_df$annotation
Yu_cluster_annotations_df$annotation <- NULL

# View the dataframe
head(Yu_cluster_annotations_df)
#                     Yu_celltype
# Human_Endo NOSTRIN  Endothelial
# Human_LAMP5 ABO             ExN
# Human_Astro_1 FGFR3   Astrocyte
# Human_PVALB ADAMTS5         InN
# Human_LAMP5 COL25A1         ExN
# Human_SOX11 EBF2            ExN




# === Yu region ===
# Replace '-' with '_' in the levels of sce.yu$space_anno
levels(sce.yu$space_anno) <- gsub("-", "_", levels(sce.yu$space_anno))
levels(sce.yu$space_anno) <- gsub("/", "_", levels(sce.yu$space_anno))
levels(sce.yu$space_anno) <- gsub("Cortical interneuron_like", "interneuron_like", levels(sce.yu$space_anno))

# Create a dataframe from the SCE object
cluster_annotations_df <- data.frame(
  annotation = sce.yu$ident,
  Yu_region = sce.yu$space_anno
)

# Remove duplicate rows to get unique mappings
Yu_region_annotations_df <- unique(cluster_annotations_df)

# Set the 'annotation' column as rownames
rownames(Yu_region_annotations_df) <- Yu_region_annotations_df$annotation
Yu_region_annotations_df$annotation <- NULL

# View the dataframe
head(Yu_region_annotations_df)
# Yu_region
# Endo NOSTRIN        non_neuron
# LAMP5 ABO                  BLA
# Astro_1 FGFR3       non_neuron
# PVALB ADAMTS5 interneuron_like
# LAMP5 COL25A1              BLA
# SOX11 EBF2                  PL

Yu_annotations <- merge(Yu_cluster_annotations_df, Yu_region_annotations_df, by="row.names")
# Set the "row.names" column as the actual row names of the merged dataframe
rownames(Yu_annotations) <- Yu_annotations$Row.names

# Remove the "row.names" column
Yu_annotations$Row.names <- NULL

# View the merged dataframe
head(Yu_annotations)
# Yu_celltype  Yu_region
# Astro_1 FGFR3   Astrocyte non_neuron
# Astro_2 FGFR3   Astrocyte non_neuron
# Astro_3 FGFR3   Astrocyte non_neuron
# Astro_4 FGFR3   Astrocyte non_neuron
# CALCR LHX8            InN    COA_MEA
# DRD2 ISL1             InN        CEA


# Define colors for each celltype
celltype_colors <- list(
  Totty_celltype = c(
    Excit = "firebrick3",
    Inhib = "dodgerblue3",
    Non_neuronal = "black"
  ),
  Yu_celltype = c(
    ExN = "firebrick3",
    InN = "dodgerblue3",
    Astrocyte = "#000000",
    Endothelial = "#444444",
    Microglia = "#666666",
    OPC = "#888888",
    Oligodendrocyte = "#999999"
  ),
  Yu_region = c(
    non_neuron = "black",
    BLA = "firebrick1",
    interneuron_like = "blue",
    IA = "darkgreen",
    COA_MEA = "gold1",
    CEA = "gold2",
    PL ="gold3",
    NLOT = "goldenrod"
  )
)


pheatmap(Totty_Yu_cor_matrix,
         color = colorRampPalette(c("#792b87", "white", "#226b38"))(50),
         main = "Pearson Correlations - Totty et al vs Yu et al",
         annotation_row = Totty_cluster_annotations_df,
         annotation_col = Yu_annotations,
         annotation_colors = celltype_colors,
         #scale="column",
         cutree_rows = 6, cutree_cols = 6,
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename=here(plot_dir, "Pearson_Totty_vs_Yu_pheatmap.pdf"),
         width=12,
         height=8)

pheatmap(Mac_Totty_cor_matrix,
         color = colorRampPalette(c("#792b87", "white", "#226b38"))(50),
         main = "Pearson Correlations - Macaque vs Human BLA",
         annotation_col = Totty_cluster_annotations_df,
         annotation_colors = celltype_colors,
         cutree_rows = 7, cutree_cols = 7,
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename=here(plot_dir, "Pearson_Mac_vs_Totty_pheatmap.pdf"),
         width=12,
         height=8)

pheatmap(Mac_Yu_cor_matrix,
         color = colorRampPalette(c("#792b87", "white", "#226b38"))(50),
         main = "Pearson Correlations - Macaque vs Human AMY (Yu et al 2023)",
         annotation_col = Yu_annotations,
         annotation_colors = celltype_colors,
         cutree_rows = 6, cutree_cols = 6,
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename=here(plot_dir, "Pearson_Mac_vs_Yu_pheatmap.pdf"),
         width=12,
         height=8)



pheatmap(Bab_Mac_cor_matrix,
         color = colorRampPalette(c("#792b87", "white", "#226b38"))(50),
         main = "Pearson Correlations - Baboon vs Macaque",
         #scale="column",
         cutree_rows = 6, cutree_cols = 6,
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename=here(plot_dir, "Pearson_Bap_vs_Mac_pheatmap.pdf"),
         width=12,
         height=8)

pheatmap(Bab_Totty_cor_matrix,
         color = colorRampPalette(c("#792b87", "white", "#226b38"))(50),
         main = "Pearson Correlations - Baboon vs Human BLA",
         annotation_col = Totty_cluster_annotations_df,
         annotation_colors = celltype_colors,
         cutree_rows = 7, cutree_cols = 7,
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename=here(plot_dir, "Pearson_Bab_vs_Totty_pheatmap.pdf"),
         width=12,
         height=8)

pheatmap(Bab_Yu_cor_matrix,
         color = colorRampPalette(c("#792b87", "white", "#226b38"))(50),
         main = "Pearson Correlations - Baboon vs Human AMY (Yu et al 2023)",
         annotation_col = Yu_annotations,
         annotation_colors = celltype_colors,
         cutree_rows = 6, cutree_cols = 6,
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename=here(plot_dir, "Pearson_Bab_vs_Yu_pheatmap.pdf"),
         width=12,
         height=8)


