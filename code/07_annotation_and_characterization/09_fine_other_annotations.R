
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
plot_dir = here("plots", "07_annotation","non-neurons")
#processed_dir = here("processed-data","07_annotation")
processed_dir <- here("processed-data")

# load sce
sce <- readRDS(here("processed-data", "JHPCE","sce_broad_annotations.rds"))
sce

colnames(colData(sce))
# [1] "Sample_num"             "Sample"                 "Species"                "Subject"               
# [5] "Sex"                    "Region"                 "Subregion"              "DV_axis"               
# [9] "PI.NeuN"                "orig.ident"             "nCount_originalexp"     "nFeature_originalexp"  
# [13] "batch"                  "Barcode"                "sum"                    "detected"              
# [17] "subsets_Mito_sum"       "subsets_Mito_detected"  "subsets_Mito_percent"   "total"                 
# [21] "high_mito"              "low_lib"                "low_genes"              "discard_auto"          
# [25] "doubletScore"           "sizeFactor"             "species"                "unintegrated_clusters" 
# [29] "seurat_clusters"        "integrated_snn_res.0.5" "ident"                  "key"                   
# [33] "broad_celltype"   


unique(sce$broad_celltype)
# [1] "Inhibitory"   "Excitatory"   "Non-neuronal"

# ======== Heatmap of canonical marker genes ========
# get only non-neurons
sce.nn <- sce[, sce$broad_celltype == "Non-neuronal"]

genes <- c("SLC1A2", # astrocyte
           "FOXJ1", # ependymal
           "MBP", # oligodendrocyte
           "TMEM119",  # microglia
           "RGS5", # Endothelial
           "CD9" # OPC
           )

# only keep features that are in the sce
genes <- genes[genes %in% rownames(sce)]

genes
plotDots(sce.nn,features=genes,group="ident", center=TRUE, scale=TRUE)


# drop ident 33, merge 3,6,14,23; merge 4 and 5
sce.nn$ident[sce.nn$ident == 33] <- NA
sce.nn$ident[sce.nn$ident %in% c(3,6,14,23)] <- 3
sce.nn$ident[sce.nn$ident %in% c(4,5)] <- 4

# drop NA and re-factor levels
sce.nn$ident <- droplevels(sce.nn$ident)
sce.nn <- sce.nn[, !is.na(sce.nn$ident)]


plotDots(sce.nn,features=genes,group="ident", center=TRUE, scale=TRUE)

# rename 3 to Oligodendrocyte, 4 to Astrocyte, 7 to Microglia, 30 to Ependymal, 31 to OPC, and 32 to Endothelial
sce.nn$ident <- as.character(sce.nn$ident)
sce.nn$ident[sce.nn$ident == 3] <- "Oligodendrocyte"
sce.nn$ident[sce.nn$ident == 4] <- "Astrocyte"
sce.nn$ident[sce.nn$ident == 7] <- "Microglia"
sce.nn$ident[sce.nn$ident == 30] <- "Ependymal"
sce.nn$ident[sce.nn$ident == 31] <- "OPC"
sce.nn$ident[sce.nn$ident == 32] <- "Endothelial"

plotReducedDim(sce.nn, dimred = "UMAP", colour_by = "ident")

sce.nn$fine_celltype <- sce.nn$ident

# save
saveRDS(sce.nn, file = here("processed-data", "JHPCE", "sce_other_final_celltypes.rds"))
