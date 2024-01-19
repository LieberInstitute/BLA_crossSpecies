# srun -n 1 --mem=124G --cpus-per-task=1 --pty bash -i
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

sce <- readRDS(here("processed-data", "07_annotation", "sce_inhib_subclustering.rds"))
colnames(colData(sce))
sce.inhib <- sce


methods <- c("louvain", "leiden", "walktrap")

ks <- c(15, 20, 25, 30, 35 , 40, 45, 50, 60)

for (i in 1:length(methods)) {
    
    # get method name
    method <- methods[i]
    
    plots <- list()
    for (j in 1:length(ks)) {
        cluster <-  paste0("k.", ks[j], "_cluster.fun.", method)
        print(cluster)
        plots[[j]] <- plotReducedDim(sce.inhib, dimred="UMAP", colour_by = cluster, point_size=0.1) +
            ggtitle(cluster) +
            theme(legend.position="none")
    }
    
    pdf(here(plot_dir, paste0(methods[i], "_subclustering.pdf")), width = 20, height = 20)
    print((plots[[1]] + plots[[2]] + plots[[3]]) / (plots[[4]] + plots[[5]] + plots[[6]]) / (plots[[7]] + plots[[8]] + plots[[9]]))
    dev.off()
    

}



pdf(here(plot_dir, "leiden_subclustering.pdf"), width = 20, height = 20)
plotReducedDim(sce.inhib, dimred="UMAP", colour_by = "k.60_cluster.fun.louvain", point_size=0.1)
dev.off()



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


# === leiden clustering ===
# excitatory cells 
excit.leiden.k75 <- clusterCells(sce.excit, 
                             use.dimred="PCA", 
                             BLUSPARAM=NNGraphParam(cluster.fun="leiden", k=75)
                             )

excit.leiden.k100 <- clusterCells(sce.excit, 
                                 use.dimred="PCA", 
                                 BLUSPARAM=NNGraphParam(cluster.fun="leiden", k=100)
)

# inhibitory cells
inhib.leiden <- clusterCells(sce.inhib, 
                             use.dimred="PCA", 
                             BLUSPARAM=NNGraphParam(cluster.fun="leiden", k=50)
                             )

# other cells
cother.leiden.k75 <- clusterCells(sce.other, 
                             use.dimred="PCA", 
                             BLUSPARAM=NNGraphParam(cluster.fun="leiden", k=75)
                             )

cother.leiden.k100 <- clusterCells(sce.other, 
                                 use.dimred="PCA", 
                                 BLUSPARAM=NNGraphParam(cluster.fun="leiden", k=100)
)

# add back to sce 
sce.excit$leiden_k75 <- excit.leiden.k75
sce.excit$leiden_k100 <- excit.leiden.k100

sce.inhib$leiden_k50 <- inhib.leiden

sce.other$leiden_k75 <- cother.leiden.k75
sce.other$leiden_k100 <- cother.leiden.k100

# === save SCE ===
#saveRDS(sce.excit, file = here(processed_dir, "sce.excit.leiden.rds"))
#saveRDS(sce.inhib, file = here(processed_dir, "sce.inhib.leiden.rds"))
#saveRDS(sce.other, file = here(processed_dir, "sce.other.leiden.rds"))
 
# =========== Run UMAP and Plot ===========

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


plotReducedDim(sce.excit, dimred = "UMAP", colour_by = "leiden_k75", point_size=0.1) +
  ggtitle("Excitatory: Leiden k=75")
plotReducedDim(sce.excit, dimred = "UMAP", colour_by = "leiden_k100", point_size=0.1) +
  ggtitle("Excitatory: Leiden k=100")

plotReducedDim(sce.inhib, dimred = "UMAP", colour_by = "leiden_k50", point_size=0.1) +
  ggtitle("Inhibitory: Leiden k=50")

plotReducedDim(sce.other, dimred = "UMAP", colour_by = "leiden_k75", point_size=0.1) + 
  ggtitle("Other: Leiden k=75")
plotReducedDim(sce.other, dimred = "UMAP", colour_by = "leiden_k100", point_size=0.1) + 
    ggtitle("Other: Leiden k=100")



# === Plot UMAPs ===

# excitatory cells
pdf(here(plot_dir, "UMAP_excitatory.pdf"), width = 15, height = 5)
p1 <- plotReducedDim(sce.excit, dimred = "UMAP", colour_by = "leiden_k75", point_size=0.1)
p2 <- plotReducedDim(sce.excit, dimred = "UMAP", colour_by = "Species", point_size=0.1)
p3 <- plotReducedDim(sce.excit, dimred = "UMAP", colour_by = "Subregion", point_size=0.1)
p1+p2+p3
dev.off()

# inhibitory cells
pdf(here(plot_dir, "UMAP_inhibitory.pdf"), width = 15, height = 5)
p1 <- plotReducedDim(sce.inhib, dimred = "UMAP", colour_by = "leiden_k50", point_size=0.1)
p2 <- plotReducedDim(sce.inhib, dimred = "UMAP", colour_by = "Species", point_size=0.1)
p3 <- plotReducedDim(sce.inhib, dimred = "UMAP", colour_by = "Subregion", point_size=0.1)
p1+p2+p3
dev.off()

# other cells
pdf(here(plot_dir, "UMAP_other.pdf"), width = 15, height = 5)
p1 <- plotReducedDim(sce.other, dimred = "UMAP", colour_by = "leiden_k100", point_size=0.1)
p2 <- plotReducedDim(sce.other, dimred = "UMAP", colour_by = "Species", point_size=0.1)
p3 <- plotReducedDim(sce.other, dimred = "UMAP", colour_by = "Subregion", point_size=0.1)
p1+p2+p3
dev.off()


# ========== Marker genes ===========
# lets try to figure out how decent these clustering results are

# dot plots for inhibitory markers

features <- c("GAD1", "GAD2",
              "SST", "PVALB",
              "LAMP5", "VIP",
              "NPY", "CORT",
              "CCK","CRH",
              "PDYN", "PENK",
              "PRKCD","CRHBP",
              "GAL", "DRD1", 
              "DRD2", "FOXP2",
              "TSHZ1", "NTS", 
              "NXPH1","CARTPT",
              "SNCG", "NOS1")

# drop features not in sce.inhib
features <- features[features %in% rownames(sce.inhib)]
features

plotDots(sce.inhib, 
         features = features, 
         group = "leiden_k50", 
         center=TRUE, 
         scale=TRUE)

# annotations
# 1 = LAMP5_NOS1
# 2 = VIP_1
# 3 = CCK_1
# 4 = NTS
# 5 = TSHZ1_1
# 6 = SST_1
# 7 = VIP_2
# 8 = PVALB_CRHBP
# 9 = PENK_DRD2
# 10 = LAMP5_2
# 11 = PVALB_2
# 12 = PVABL_DRD2
# 13 = SST_2
# 14 = VIP_3
# 15 = TSHZ1_2
# 16 = TSHZ1_PRKCD

# make new $fine_type annotation column 
sce.inhib$fine_type <- c()

# assign annotations
sce.inhib$fine_type[sce.inhib$leiden_k50 == 1] <- "LAMP5_NOS1"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 2] <- "VIP_1"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 3] <- "CCK_1"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 4] <- "NTS"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 5] <- "TSHZ1_1"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 6] <- "SST_1"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 7] <- "VIP_2"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 8] <- "PVALB_CRHBP"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 9] <- "PENK_DRD2"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 10] <- "LAMP5_2"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 11] <- "PVALB_2"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 12] <- "PVABL_DRD2"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 13] <- "SST_2"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 14] <- "VIP_3"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 15] <- "TSHZ1_2"
sce.inhib$fine_type[sce.inhib$leiden_k50 == 16] <- "TSHZ1_PRKCD"

sce.inhib$fine_type <- as.factor(sce.inhib$fine_type)

unique(sce.inhib$fine_type)

# remake umap with annotations
pdf(here(plot_dir, "UMAP_inhib_annotated.pdf"), width = 5, height = 5)
plotReducedDim(sce.inhib, dimred="UMAP", colour_by = "fine_type", point_size=0.1, text_by="fine_type")
dev.off()


# other cells
pdf(here(plot_dir, "UMAP_inhib_annotated_markers.pdf"), width = 15, height = 15)
features2plots <- c("SST", "NOS1", "PVALB", "LAMP5", "VIP", "PENK", "TSHZ1", "NTS", "PRKCD")
# intialize variable to hold plots in loop
plots <- list()
for (i in features2plots) {
    plots[[i]] <- plotReducedDim(sce.inhib, dimred="UMAP", colour_by = i, point_size=0.1) +
        scale_color_gradient(low = "grey", high = "red") +
        ggtitle(i)
}

# plot grid of 9 plots
gridExtra::grid.arrange(grobs = plots, ncol = 3)
dev.off()




# ===== testing


inhib.out <- clusterSweep(reducedDim(sce.inhib, "PCA"), 
                          NNGraphParam(), 
                          k=as.integer(c(10, 50, 100)),
                          cluster.fun=c("leiden"),
                          BPPARAM=BiocParallel::MulticoreParam(8))

colData(sce.inhib) <- cbind(colData(sce.inhib), inhib.out$clusters)

