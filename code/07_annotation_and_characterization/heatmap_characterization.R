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
library(MetaNeighbor)
library(tidyverse)
library(ggplot2)
library(ghibli)
library(RColorBrewer)
library(ComplexHeatmap)

# directories
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "08_species_comparisons")

sce.inhib <- readRDS(here(processed_dir, "sce.inhib.final.rds"))
sce.excit <- readRDS(here(processed_dir, "sce.excit.integrated.annotated.rds"))


# join SCE objects
sce <- cbind(sce.inhib, sce.excit)
sce
# class: SingleCellExperiment 
# dim: 13874 99588 
# metadata(0):
#     assays(2): counts logcounts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(0):
#     colnames(99588): AAACCCAAGCTAAATG-1 AAACCCAGTTATCCAG-1 ... TGGGAAGGTTAGCGGA-1
# TTCACGCGTAGTCTTG-1
# colData names(61): orig.ident nCount_originalexp ... k.60_cluster.fun.louvain fine_type
# reducedDimNames(2): PCA UMAP
# mainExpName: originalexp
# altExpNames(0):

colnames(colData(sce))
# [1] "orig.ident"                "nCount_originalexp"        "nFeature_originalexp"     
# [4] "Sample_num"                "Sample"                    "Species"                  
# [7] "Subject"                   "Sex"                       "Region"                   
# [10] "Subregion"                 "DV_axis"                   "PI.NeuN"                  
# [13] "batch"                     "Barcode"                   "sum"                      
# [16] "detected"                  "subsets_Mito_sum"          "subsets_Mito_detected"    
# [19] "subsets_Mito_percent"      "total"                     "high_mito"                
# [22] "low_lib"                   "low_genes"                 "discard_auto"             
# [25] "doubletScore"              "sizeFactor"                "species"                  
# [28] "originalexp_snn_res.2"     "seurat_clusters"           "integrated_snn_res.0.5"   
# [31] "ident"                     "key"                       "broad_celltype"           
# [34] "k.15_cluster.fun.walktrap" "k.20_cluster.fun.walktrap" "k.25_cluster.fun.walktrap"
# [37] "k.30_cluster.fun.walktrap" "k.35_cluster.fun.walktrap" "k.40_cluster.fun.walktrap"
# [40] "k.45_cluster.fun.walktrap" "k.50_cluster.fun.walktrap" "k.60_cluster.fun.walktrap"
# [43] "k.15_cluster.fun.leiden"   "k.20_cluster.fun.leiden"   "k.25_cluster.fun.leiden"  
# [46] "k.30_cluster.fun.leiden"   "k.35_cluster.fun.leiden"   "k.40_cluster.fun.leiden"  
# [49] "k.45_cluster.fun.leiden"   "k.50_cluster.fun.leiden"   "k.60_cluster.fun.leiden"  
# [52] "k.15_cluster.fun.louvain"  "k.20_cluster.fun.louvain"  "k.25_cluster.fun.louvain" 
# [55] "k.30_cluster.fun.louvain"  "k.35_cluster.fun.louvain"  "k.40_cluster.fun.louvain" 
# [58] "k.45_cluster.fun.louvain"  "k.50_cluster.fun.louvain"  "k.60_cluster.fun.louvain" 
# [61] "fine_type"    


# normalize
factors <- calculateSumFactors(sce, cluster = sce$fine_type)
sce <- logNormCounts(sce, size.factors=factors)

# get HVGs
model <- modelGeneVar(sce)
hvg <- getTopHVGs(model, n=3000)

# 2. Extract Expression Data for HVGs
out <- aggregateAcrossCells(sce,
                            ids=sce$fine_type, 
                            subset.row = hvg,
                            use.assay.type = "logcounts"
)

distance_matrix <- dist(t(logcounts(out)))
hc <- hclust(distance_matrix, method = "ward.D2")
plot(hc)

# 4. Save the Column Dendrogram
column_dendrogram <- as.dendrogram(hc)
plot(column_dendrogram)


heatmap_annotation <- HeatmapAnnotation(df = data.frame(sce$total),
                                    which = "column") 
                                                                         
