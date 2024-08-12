
library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(DeconvoBuddies)
library(dplyr)
library(bluster)
library(pheatmap)
library(ggpubr)


## save directories
plot_dir = here("plots", "07_annotation","05_inhib_annotations")
#processed_dir = here("processed-data","07_annotation")
processed_dir <- here("processed-data")

# load sce
sce <- readRDS(here("processed-data", "JHPCE","sce_inhib_final_subclusters.rds"))
sce


# ==== New subclustering =====

# inhibatory clustering
set.seed(112)
inhib.out <- clusterCells(sce, use.dimred="PCA", 
                          BLUSPARAM=NNGraphParam(cluster.fun="leiden", k=50), full=TRUE)

sce$new_subclusters <- inhib.out$clusters

png(here(plot_dir,"UMAP_inhib_subclusters_k50.png"), width=10, height=10, units="in", res=300)
plotReducedDim(sce, dimred="UMAP", colour_by="new_subclusters", point_size=0.75)
dev.off()




# ====== Marker genes ======

markers <- findMarkers_1vAll(
    sce,
    assay_name = "logcounts",
    cellType_col = "new_subclusters",
    mod = "~species",
    verbose = TRUE
)


# Group by cellType.target and filter top 2 genes
top_genes <- markers %>%
    group_by(cellType.target) %>%
    top_n(n = 2, wt = -rank_marker) %>%
    arrange(cellType.target, rank_marker) %>%
    summarise(gene1_gene2 = paste(gene[1], gene[2], sep = "_"))

# Display the result
print(top_genes)

# named vector for mapping
#name_map <- setNames(lookup_table$gene1_gene2, lookup_table$cellType.target)

#rename fine celltypes
#sce$fine_celltype <- name_map[as.character(sce$fine_celltype)]
#unique(sce$fine_celltype)
# [1] "ESR1_ADRA1A"     "MEIS2_COL25A1"   "MEIS1_PARD3B"    "ZBTB20_SLC4A4"   "PEX5L_MYRIP"     "ST18_ABCA8"      "ADARB2_TRPS1"   
# [8] "GULP1_TRHDE"     "RXFP1_KIAA1217"  "GRIK3_TNS3"      "SATB2_MPPED1"    "SLC17A8_ST8SIA2"


# ======== Heatmap of canonical marker genes ========

genes <- c("SST", "VIP","PVALB","CCK","LAMP5","CARTPT","FOXP2","TSHZ1","CRH", "NOS1","PPP1R1B","DRD1","DRD2","PRKCD")

plotDots(sce,features=genes,group="new_subclusters", center=TRUE, scale=TRUE)

subcluster_names <- c("LAMP5_NTNG1", "VIP_ADRA1B", "THSD7B_CALB2", "TSHZ1.1",
                      "SST_NOS1", "PVALB_MYO5B", "PRKCD_DRD2", "SST_PRKCQ",
                      "LAMP5_KIT", "CCK_CNR1", "ZFHX3_SCN5A", "PVALB_UNC5B",
                      "CARTPT_CDH23", "ST18_IL1RAPL2", "LHX8_ANGPT1", "VIP_PLPP4",
                      "PRKCD_DRD1", "TSHZ1.2")

# Create a named vector for mapping
names(subcluster_names) <- 1:length(subcluster_names)

sce$new_subclusters <- as.character(sce$new_subclusters)

# Update new_subclusters using the named vector
sce$new_subclusters <- subcluster_names[as.character(sce$new_subclusters)]

# Verify the updates
unique(sce$new_subclusters)

sce$fine_celltype <- sce$new_subclusters


# ===== save annotated subclusters ======#

# clean up SCE

# drop k.25_cluster.fun.leiden
colData(sce)["inhib.out$clusters"] <- NULL
sce$new_subclusters <- NULL

saveRDS(sce, here("processed-data", "sce_inhib_final_subclusters_annotated.rds"))




# ====== Heatmap of canonical marker genes ======

genes <- c("CARTPT","CCK", "LAMP5","PRKCD", "DRD1","DRD2", "PVALB","SST","NOS1",
           "FOXP2", "VIP","CRH")

plotGroupedHeatmap(sce, features=genes, 
                   group="fine_celltype", 
                   center=TRUE, 
                   scale=TRUE,
                   cluster_rows=FALSE,
                   cluster_cols=FALSE
                   ) 
    

    
# ======= ComplexHeatmap ========
library(ComplexHeatmap)


# Extract the expression data for the specified genes and fine cell types
expr_data <- logcounts(sce)[genes, ]
cell_types <- sce$fine_celltype

# Optionally, center and scale the data
expr_data_scaled <- t(scale(t(expr_data), center = TRUE, scale = TRUE))

# Prepare annotations (if any)
ha <- HeatmapAnnotation(
    df = data.frame(fine_celltype = cell_types),
    col = list(fine_celltype = c("1" = "red", "2" = "blue", "3" = "green", "4" = "purple",
                                 "5" = "orange", "6" = "brown", "7" = "pink", "8" = "cyan",
                                 "9" = "yellow", "10" = "black", "11" = "grey", "12" = "lightblue",
                                 "13" = "lightgreen", "14" = "lightcoral", "15" = "lightgoldenrod",
                                 "16" = "lightpink", "17" = "lightsalmon", "18" = "lightseagreen"))
)

# Create the heatmap
heatmap <- Heatmap(expr_data_scaled,
                   name = "Expression",
                   row_names_side = "left",
                   column_names_side = "bottom",
                   top_annotation = ha,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   show_row_names = TRUE,
                   show_column_names = FALSE)

# Draw the heatmap
draw(heatmap)
