
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
plot_dir = here("plots", "07_annotation","excit_subclusters")
#processed_dir = here("processed-data","07_annotation")
processed_dir <- here("processed-data")

# load sce
sce <- readRDS(here("processed-data", "JHPCE","sce_excit_final_subclusters.rds"))
sce
# class: SingleCellExperiment 
# dim: 13842 76519 
# metadata(0):
# assays(2): counts logcounts
# rownames(13842): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(0):
# colnames(76519): AAACCCACAGGTCCCA-1 AAACCCATCCATTGTT-1 ... TTCACGCGTAGTCTTG-1 TTCCGGTCACATCCCT-1
# colData names(35): orig.ident nCount_originalexp ... integrated_snn_res.0.4 subcluster_idents
# reducedDimNames(2): PCA UMAP
# mainExpName: NULL
# altExpNames(0):


# ==== New subclustering =====

methods <- c("leiden")
ks <- c(100)

# Excitatory clustering
excit.out <- clusterCells(sce, use.dimred="PCA", 
             BLUSPARAM=NNGraphParam(cluster.fun="leiden", k=100), full=TRUE)

colData(sce) <- cbind(colData(sce), excit.out$clusters)

png(here(plot_dir,"UMAP_excit_subclusters_leiden_k100.png"), width=10, height=10, units="in", res=300)
plotReducedDim(sce, dimred="UMAP", colour_by="k.100_cluster.fun.leiden")
dev.off()



g <- excit.out$objects$graph
ratio <- pairwiseModularity(g, sce$k.100_cluster.fun.leiden, as.ratio=TRUE)
dim(ratio)


png(here(plot_dir,"Heatmap_pairwise_similarity_leiden_k100.png"), width=10, height=10, units="in", res=300)
p1 <- pheatmap(log2(ratio+1), cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100))
p1
dev.off()

# combine cluster 13 -> 7
sce$k.100_cluster.fun.leiden[sce$k.100_cluster.fun.leiden == 13] <- 7

png(here(plot_dir,"UMAP_excit_subclusters_leiden_k100_dropped_cluster13.png"), width=10, height=10, units="in", res=300)
plotReducedDim(sce, dimred="UMAP", colour_by="k.100_cluster.fun.leiden")
dev.off()


# ====== Marker genes ======

markers <- findMarkers_1vAll(
    sce,
    assay_name = "logcounts",
    cellType_col = "k.100_cluster.fun.leiden",
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
name_map <- setNames(lookup_table$gene1_gene2, lookup_table$cellType.target)

#rename fine celltypes
sce$fine_celltype <- name_map[as.character(sce$fine_celltype)]
unique(sce$fine_celltype)
# [1] "ESR1_ADRA1A"     "MEIS2_COL25A1"   "MEIS1_PARD3B"    "ZBTB20_SLC4A4"   "PEX5L_MYRIP"     "ST18_ABCA8"      "ADARB2_TRPS1"   
# [8] "GULP1_TRHDE"     "RXFP1_KIAA1217"  "GRIK3_TNS3"      "SATB2_MPPED1"    "SLC17A8_ST8SIA2"


# ===== save annotated subclusters ======#

# clean up SCE

# drop k.25_cluster.fun.leiden
sce$k.25_cluster.fun.leiden <- NULL
sce$k.50_cluster.fun.leiden <- NULL
sce$k.75_cluster.fun.leiden <- NULL
sce$k.100_cluster.fun.leiden <- NULL
sce$subcluster_idents <- NULL

saveRDS(sce, here("processed-data", "sce_excit_final_subclusters_annotated.rds"))


