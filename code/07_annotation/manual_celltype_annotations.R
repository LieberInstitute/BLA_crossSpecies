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

## save directories
plot_dir = here("plots", "07_annotation", "cluster_check")
processed_dir = here("processed-data","07_annotation", "cluster_check")

# load sce
# save combined, uncorrected sce
load(here("processed-data","06_clustering", "seurat_integrated_final.rda"))
seurat.int

Assays(seurat.int)
# [1] "originalexp" "integrated" 

# convert to sce
sce <- as.SingleCellExperiment(seurat.int, assay = "originalexp")
sce
# class: SingleCellExperiment 
# dim: 13874 177346 
# metadata(0):
#     assays(2): counts logcounts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(0):
#     colnames(177346): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 35_TTTGGTTGTGGCTAGA-1 35_TTTGTTGCACATGGTT-1
# colData names(23): orig.ident nCount_originalexp ...
# integrated_snn_res.0.5 ident
# reducedDimNames(2): PCA UMAP
# mainExpName: originalexp
# altExpNames(0):


unique(colnames(colData(sce)))
# [1] "orig.ident"             "nCount_originalexp"     "nFeature_originalexp"  
# [4] "batch"                  "Sample"                 "Barcode"               
# [7] "sum"                    "detected"               "subsets_Mito_sum"      
# [10] "subsets_Mito_detected"  "subsets_Mito_percent"   "total"                 
# [13] "high_mito"              "low_lib"                "low_genes"             
# [16] "discard_auto"           "doubletScore"           "sizeFactor"            
# [19] "species"                "originalexp_snn_res.2"  "seurat_clusters"       
# [22] "integrated_snn_res.0.5" "ident"   




# ========= Checking for clusters of doublets ========

pdf(here(plot_dir, "UMAP_doublet_scores.pdf"))
plotReducedDim(sce, dimred = "UMAP", colour_by = "doubletScore")
dev.off()

pdf(here(plot_dir, "UMAP_clusters.pdf"))
plotReducedDim(sce, dimred = "UMAP", colour_by = "ident")
dev.off()


# ========== plotting UMAPs of broad marker genes ========

# plot some broad marker genes
pdf(here(plot_dir, "UMAP_broad_markers.pdf"), width=25, height=10)
p1 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "SNAP25")
p2 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "SLC17A7")
p3 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "GAD2")
p4 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "MBP")
p5 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "GFAP")
p6 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "ident")

print((p1 | p2 | p3) / (p4 | p5 | p6))
dev.off()



# =========== plotting heatmap of broad marker genes =========

features <- c("SNAP25", "SYT1",
              "SLC17A7","SLC17A6",
              "GAD1", "GAD2",
              "LAMP5", "VIP",
              "FOXP2", "MEIS2",
              "SST", "PVALB",
              "MBP", "MOBP",
              "GFAP", "AQP4",
              "CD74","C3",
              "PDGFRA", "VCAN",
              "CLDN5", "FLT1"
            )

# drop any genes in features that are not in sce
features <- features[features %in% rownames(sce)]
features
# [1] "SNAP25"  "SLC17A6" "GAD1"    "GAD2"    "LAMP5"   "VIP"     "FOXP2"  
# [8] "MEIS2"   "SST"     "PVALB"   "MBP"     "GFAP"    "CD74"    "VCAN"   
# [15] "CLDN5"   "FLT1"  

pdf(here(plot_dir, "Heatmap_broad_markers.pdf"), width=10, height=10)
plotGroupedHeatmap(sce, 
                   group = "ident",
                   features = features,
                   symmetric=TRUE,
                   center=TRUE
                   )
dev.off()

# dot plots
pdf(here(plot_dir, "Dotplot_broad_markers.pdf"), width=10, height=10)
plotDots(sce, 
         group = "ident",
         features = features
        )
dev.off()


# ======== Annotation broad cell types ========

excit <- c(27,14,19,33,1,13,16,7,10,2,9,4,25)
inhib <- c(17,15,21,18,8,24, 11 ,12, 21, 34, 20)
non_neuronal <- c(0, 3,5,6,22,23,26,28,29,30,31,32)

sce$ident <- as.factor(sce$ident)

# create new broad_celltype column based on above clusters
colData(sce)$broad_celltype <- "NA"
colData(sce)$broad_celltype[which(colData(sce)$ident %in% excit)] <- "Excitatory"
colData(sce)$broad_celltype[which(colData(sce)$ident %in% inhib)] <- "Inhibitory"
colData(sce)$broad_celltype[which(colData(sce)$ident %in% non_neuronal)] <- "Non-neuronal"

# plot UMAP with broad cell types
pdf(here(plot_dir, "UMAP_broad_celltype_labels.pdf"))
plotReducedDim(sce, dimred = "UMAP", colour_by = "broad_celltype")
dev.off()


# ========== Comparing macaque, baboon, and human cells ========

#subset to different species
sce_macaque <- sce[,which(colData(sce)$species == "macaque")]
sce_baboon <- sce[,which(colData(sce)$species == "baboon")]
sce_human <- sce[,which(colData(sce)$species == "human")]

pdf(here(plot_dir, "UMAP_across_species.pdf"), width=30, height=10)

p1 <- plotReducedDim(sce, dimred = "UMAP", colour_by = "broad_celltype") +
    theme_void() +
    theme(legend.position="none",
          plot.title = element_text("All Species"))

p2 <- plotReducedDim(sce_macaque, dimred = "UMAP", colour_by = "broad_celltype") +
    theme_void() +
    theme(legend.position="none",
          plot.title = element_text("Macaque"))

p3 <- plotReducedDim(sce_baboon, dimred = "UMAP", colour_by = "broad_celltype") +
    theme_void() +
    theme(legend.position="none",
          plot.title = element_text("Baboon"))

p4 <- plotReducedDim(sce_human, dimred = "UMAP", colour_by = "broad_celltype") +
    theme_void() +
    theme(legend.position="none",
          plot.title = element_text("Human"))


print((p1 | p2 | p3 | p4))
dev.off()