library("dplyr")
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
library("pheatmap")
library(tidySingleCellExperiment)

# directories
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "01_pearson_correlations")

sce.inhib <- readRDS(here(processed_dir, "sce.inhib.final.rds"))
sce.excit <- readRDS(here(processed_dir, "sce.excit.integrated.annotated.rds"))

sce <- cbind(sce.inhib, sce.excit)
sce
# class: SingleCellExperiment 
# dim: 13874 99588 
# metadata(0):
#     assays(2): counts logcounts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(0):
#     colnames(99588): AAACCCAAGCTAAATG-1 AAACCCAGTTATCCAG-1 ... TGGGAAGGTTAGCGGA-1 TTCACGCGTAGTCTTG-1
# colData names(61): orig.ident nCount_originalexp ... k.60_cluster.fun.louvain fine_type
# reducedDimNames(2): PCA UMAP
# mainExpName: originalexp
# altExpNames(0):



# ======= Subset to each species =======

sce.human <- sce[, sce$species == "human"]
sce.macaque <- sce[, sce$species == "macaque"]
sce.baboon <- sce[, sce$species == "baboon"]

# removed unused SCEs for memory
rm(sce.excit)
rm(sce.inhib)
rm(sce)

# ======= Pearson correlation within Macaque celltypes using all genes =======

library(reshape2)


# drop level ZNF804B_TTN: not found in human dataset
sce.macaque <- sce.macaque[,sce.macaque$fine_type != "ZNF804B_TTN"]
sce.human <- sce.human[,sce.human$fine_type != "ZNF804B_TTN"]
sce.baboon <- sce.baboon[,sce.baboon$fine_type != "ZNF804B_TTN"]

# drop levels
sce.macaque$fine_type <- droplevels(sce.macaque$fine_type)
sce.human$fine_type <- droplevels(sce.human$fine_type)
sce.baboon$fine_type <- droplevels(sce.baboon$fine_type)

# 1. Extract expression data
Human_exprs_data <- logcounts(sce.human)  # or use `assay(sce, "normalized")` or similar depending on your data
Mac_exprs_data <- logcounts(sce.macaque)
Bab_exprs_data <- logcounts(sce.baboon)

# 3. Subset each dataset
#Yu_exprs_data_common <- Yu_exprs_data[common_genes, ]

# scale if you want to
#Mac_exprs_data <- scale(Mac_exprs_data)
#Human_exprs_data <- scale(Human_exprs_data)
#Yu_exprs_data <- scale(Yu_exprs_data)
#Bab_exprs_data <- scale(Bab_exprs_data)

# 4. Group by clusters and compute mean expression
Mac_clusters <- unique(sce.macaque$fine_type)
Human_clusters <- unique(sce.human$fine_type)
Bab_clusters <- unique(sce.baboon$fine_type)
#Yu_clusters <- unique(colData(sce.yu)$ident)


#  function to get mean expression for each cluster
compute_mean_exprs <- function(exprs_data, sce, cluster_col, cluster) {
    subset_data <- exprs_data[, colData(sce)[[cluster_col]] == cluster]
    print(paste0(cluster, ": ", dim(subset_data)))
    return(rowMeans(subset_data))
}

# Now use the function for each dataset:

# For Mac data
Mac_mean_exprs_list <- lapply(common_clusters, function(cluster) {
    compute_mean_exprs(Mac_exprs_data, sce.macaque, "fine_type", cluster)
})

# For Human data
Human_mean_exprs_list <- lapply(Human_clusters, function(cluster) {
    compute_mean_exprs(Human_exprs_data, sce.human, "fine_type", cluster)
})

# For Yu data
# Yu_mean_exprs_list <- lapply(Yu_clusters, function(cluster) {
#     compute_mean_exprs(Yu_exprs_data_common, sce.yu, "ident", cluster)
# })

# For Baboon data
Bab_mean_exprs_list <- lapply(Bab_clusters, function(cluster) {
    compute_mean_exprs(Bab_exprs_data, sce.baboon, "fine_type", cluster)
})


# Convert the list of mean expressions to a matrix
Mac_mean_exprs_matrix <- do.call(cbind, Mac_mean_exprs_list)
colnames(Mac_mean_exprs_matrix) <- Mac_clusters

Human_mean_exprs_matrix <- do.call(cbind, Human_mean_exprs_list)
colnames(Human_mean_exprs_matrix) <- Human_clusters

#Yu_mean_exprs_matrix <- do.call(cbind, Yu_mean_exprs_list)
#colnames(Yu_mean_exprs_matrix) <- Yu_clusters

Bab_mean_exprs_matrix <- do.call(cbind, Bab_mean_exprs_list)
colnames(Bab_mean_exprs_matrix) <- Bab_clusters




# 5. Calculate Pearson correlation
Mac_Human_cor_matrix <- cor(Mac_mean_exprs_matrix, Human_mean_exprs_matrix)
Mac_Yu_cor_matrix <- cor(Mac_mean_exprs_matrix, Yu_mean_exprs_matrix)
Human_Yu_cor_matrix <- cor(Human_mean_exprs_matrix, Yu_mean_exprs_matrix)

Bab_Mac_cor_matrix <- cor(Bab_mean_exprs_matrix, Mac_mean_exprs_matrix)
Bab_Human_cor_matrix <- cor(Bab_mean_exprs_matrix, Human_mean_exprs_matrix)
Bab_Yu_cor_matrix <- cor(Bab_mean_exprs_matrix, Yu_mean_exprs_matrix)


# 5. Visualize correlation


# ========== Get annotations ==========

# === Human celltypes ===
sce.Human$celltype <- factor(sce.Human$celltype)
levels(sce.Human$celltype) <- gsub("-", "_", levels(sce.Human$celltype))

# Create a dataframe from the SCE object
cluster_annotations_df <- data.frame(
    annotation = sce.Human$annotation,
    Human_celltype = sce.Human$celltype
)

# Remove duplicate rows to get unique mappings
Human_cluster_annotations_df <- unique(cluster_annotations_df)

# Set the 'annotation' column as rownames
rownames(Human_cluster_annotations_df) <- Human_cluster_annotations_df$annotation
Human_cluster_annotations_df$annotation <- NULL

# View the dataframe
head(Human_cluster_annotations_df)
# Human_celltype
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
    Human_celltype = c(
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


pheatmap(Human_Yu_cor_matrix,
         color = colorRampPalette(c("#792b87", "white", "#226b38"))(50),
         main = "Pearson Correlations - Human et al vs Yu et al",
         annotation_row = Human_cluster_annotations_df,
         annotation_col = Yu_annotations,
         annotation_colors = celltype_colors,
         #scale="column",
         cutree_rows = 6, cutree_cols = 6,
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename=here(plot_dir, "Pearson_Human_vs_Yu_pheatmap.pdf"),
         width=12,
         height=8)

pheatmap(Mac_Human_cor_matrix,
         color = colorRampPalette(c("#792b87", "white", "#226b38"))(50),
         main = "Pearson Correlations - Macaque vs Human BLA",
         annotation_col = Human_cluster_annotations_df,
         annotation_colors = celltype_colors,
         cutree_rows = 7, cutree_cols = 7,
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename=here(plot_dir, "Pearson_Mac_vs_Human_pheatmap.pdf"),
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

pheatmap(Bab_Human_cor_matrix,
         color = colorRampPalette(c("#792b87", "white", "#226b38"))(50),
         main = "Pearson Correlations - Baboon vs Human BLA",
         annotation_col = Human_cluster_annotations_df,
         annotation_colors = celltype_colors,
         cutree_rows = 7, cutree_cols = 7,
         show_rownames = TRUE,
         show_colnames = TRUE,
         filename=here(plot_dir, "Pearson_Bab_vs_Human_pheatmap.pdf"),
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


