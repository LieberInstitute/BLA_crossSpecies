# srun -n 1 --mem=200G --cpus-per-task=8 --pty bash -i
library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(batchelor)
library(harmony)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
library(BiocParallel)
library(bluster)

## save directories
plot_dir = here("plots", "07_annotation", "cluster_check")
processed_dir = here("processed-data","07_annotation", "cluster_check")

# save sce
sce <- readRDS(here("processed-data", "07_annotation", "sce_broad_annotations.rds"))
sce
# class: SingleCellExperiment 
# dim: 13874 177346 
# metadata(0):
#     assays(2): counts logcounts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(0):
#     colnames(177346): AAACCCAAGCTAAATG-1 AAACCCACAGGTCCCA-1 ...
# TTTGGTTGTGGCTAGA-1 TTTGTTGCACATGGTT-1
# colData names(32): Sample_num Sample ... ident key
# reducedDimNames(2): PCA UMAP
# mainExpName: originalexp
# altExpNames(0):

unique(colnames(colData(sce)))
# [1] "Sample_num"             "Sample"                 "Species"               
# [4] "Subject"                "Sex"                    "Region"                
# [7] "Subregion"              "DV_axis"                "PI.NeuN"               
# [10] "orig.ident"             "nCount_originalexp"     "nFeature_originalexp"  
# [13] "batch"                  "Barcode"                "sum"                   
# [16] "detected"               "subsets_Mito_sum"       "subsets_Mito_detected" 
# [19] "subsets_Mito_percent"   "total"                  "high_mito"             
# [22] "low_lib"                "low_genes"              "discard_auto"          
# [25] "doubletScore"           "sizeFactor"             "species"               
# [28] "originalexp_snn_res.2"  "seurat_clusters"        "integrated_snn_res.0.5"
# [31] "ident"                  "key"  

# ========= subsetting to broad celltypes ============

# subset to broad celltypes in $broad_celltype
sce.excit <- sce[, sce$broad_celltype == "Excitatory"]
sce.inhib <- sce[, sce$broad_celltype == "Inhibitory"]
sce.other <- sce[, sce$broad_celltype == "Non-neuronal"]


# ========= walktrap clustering with various k ==========

methods <- c("walktrap", "leiden", "louvain")
ks <- c(15, 20, 25, 30, 35, 40, 45, 50, 60)

# Excitatory clustering
excit.out <- clusterSweep(reducedDim(sce.excit, "PCA"), 
                    NNGraphParam(), 
                    k=as.integer(ks),
                    cluster.fun=methods,
                    BPPARAM=BiocParallel::MulticoreParam(8))

colData(sce.excit) <- cbind(colData(sce.excit), excit.out$clusters)


# inhibitory clustering
inhib.out <- clusterSweep(reducedDim(sce.inhib, "PCA"), 
                          NNGraphParam(), 
                          k=as.integer(ks),
                          cluster.fun=methods,
                          BPPARAM=BiocParallel::MulticoreParam(8))

colData(sce.inhib) <- cbind(colData(sce.inhib), inhib.out$clusters)


# non-neuronal clustering
other.out <- clusterSweep(reducedDim(sce.other, "PCA"), 
                          NNGraphParam(), 
                          k=as.integer(ks),
                          cluster.fun=methods,
                          BPPARAM=BiocParallel::MulticoreParam(8))

colData(sce.other) <- cbind(colData(sce.other), other.out$clusters)


# ========= save sce ===========
saveRDS(sce.excit, file = here("processed-data", "07_annotation", "sce_excit_subclustering.rds"))
saveRDS(sce.inhib, file = here("processed-data", "07_annotation", "sce_inhib_subclustering.rds"))
saveRDS(sce.other, file = here("processed-data", "07_annotation", "sce_other_subclustering.rds"))



# =============== Plotting ==================

# === Run UMAP ===
set.seed(111)

# excitatory cells
sce.excit <- runUMAP(sce.excit, 
                     dimred = "PCA",
                     min_dist = 0.3
)

# inhibitory cells
sce.inhib <- runUMAP(sce.inhib, 
                     dimred = "PCA",
                     min_dist = 0.3
)

# other cells
sce.other <- runUMAP(sce.other, 
                     dimred = "PCA",
                     min_dist = 0.3
)

# === Loop through UMAP plots ===

# plot UMAPs with various ks for each method. One plot per method.
# the column names follow the scheme "k.10_cluster.fun.leiden"

for (i in 1:length(methods)) {
    
    # get method name
    method <- methods[i]
    
    # initialize pdf
    pdf(here(plot_dir, paste0(methods, "_subclustering.pdf")), width = 20, height = 20)
        
    plots <- list()
    for (j in length(ks)) {
      cluster <-  paste0("k.", ks[j], "_cluster.fun.", method)
      plots[[j]] <- plotReducedDim(sce.inhib, dimred="UMAP", colour_by = cluster, point_size=0.1) +
          ggtitle(cluster)
    }
    
    gridExtra::grid.arrange(grobs = plots, ncol = 3)
    dev.off()
}








