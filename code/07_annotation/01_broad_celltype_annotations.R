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

## save directories
plot_dir = here("plots", "07_annotation", "01_broad_annotations")
processed_dir = here("processed-data","07_annotation",  "01_broad_annotations")

# load sce
load(here("processed-data","06_clustering", "seurat_integrated_final.rda"))
seurat.int

Assays(seurat.int)
# [1] "originalexp" "integrated" 

# convert to sce
sce <- as.SingleCellExperiment(seurat.int, assay = "originalexp")
sce
# class: SingleCellExperiment 
# dim: 13842 176107 
# metadata(0):
# assays(2): counts logcounts
# rownames(13842): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(0):
# colnames(176107): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
#   35_TTTGGTTGTGGCTAGA-1 35_TTTGTTGCACATGGTT-1
# colData names(23): orig.ident nCount_originalexp ...
#   integrated_snn_res.0.5 ident
# reducedDimNames(2): PCA UMAP
# mainExpName: originalexp
# altExpNames(0):


# load metadata to re-marge



unique(colnames(colData(sce)))
# [1] "orig.ident"             "nCount_originalexp"     "nFeature_originalexp"  
# [4] "batch"                  "Sample"                 "Barcode"               
# [7] "sum"                    "detected"               "subsets_Mito_sum"      
# [10] "subsets_Mito_detected"  "subsets_Mito_percent"   "total"                 
# [13] "high_mito"              "low_lib"                "low_genes"             
# [16] "discard_auto"           "doubletScore"           "sizeFactor"            
# [19] "species"                "originalexp_snn_res.2"  "seurat_clusters"       
# [22] "integrated_snn_res.0.5" "ident"   

# drop any duplicate coldata
sce <- sce[,!duplicated(colnames(colData(sce)))]

# ========= Add metadata =========

# load metadata csv
metadata <- read.csv(here("raw-data","sampleinfo", "master_sampleinfo_2023-10-09.csv"))

# drop samples with Region == DLPFC
metadata <- metadata[metadata$Region == "Amygdala",]

unique(sce$Sample)
# [1] "Br2327-3c-AMYBLA"                                    
# [2] "Br8692-4c-AMYBLA"                                    
# [3] "Br9021-5c-AMYBLA"                                    
# [4] "Br8331-34ac_scp"                                     
# [5] "Br5273-35ac_scp"                                     
# [6] "SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1"         
# [7] "SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1"          
# [8] "SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1"          
# [9] "SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1"          
# [10] "SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1"
# [11] "SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1"          
# [12] "SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1"
# [13] "VC_snRNAseq-7_LateralVentralAP1_-1"                  
# [14] "VC_snRNAseq-7_LateralVentralAP2_-2"                  
# [15] "VC_snRNAseq-7_BasalDorsalAP1_-3A"                    
# [16] "VC_snRNAseq-7_BasalDorsalAP1_-3B"                    
# [17] "VC_snRNAseq-7_BasalVentralAP1_-4"                    
# [18] "VC_snRNAseq-7_BasalVentralAP2_-5"                    
# [19] "VC_snRNAseq-7_AccBasalAP1AP2_-7"                     
# [20] "VC_snRNAseq-7_LateralDorsalAP1AP2_-8"                
# [21] "VC_snRNAseq_9_Animal4_LV2"                           
# [22] "VC_snRNAseq_9_Animal4_LV1"                           
# [23] "VC_snRNAseq_9_Animal4_LD"                            
# [24] "VC_snRNAseq_9_Animal4_L-Comb"                        
# [25] "VC_snRNAseq_9_Animal4_BV"                            
# [26] "VC_snRNAseq_9_Animal4_BD"                            
# [27] "VC_snRNAseq_9_Animal4_B-Comb"                        
# [28] "VC_snRNAseq_9_Animal4_AB"                            
# [29] "VC_snRNAseq_8_Animal3_LV"                            
# [30] "VC_snRNAseq_8_Animal3_LV3"                           
# [31] "VC_snRNAseq_8_Animal3_LD"                            
# [32] "VC_snRNAseq_8_Animal3_BV2"                           
# [33] "VC_snRNAseq_8_Animal3_BV1"                           
# [34] "VC_snRNAseq_8_Animal3_BD"                            
# [35] "VC_snRNAseq_8_Animal3_AB__"                          
# [36] "VC_snRNAseq_8_Animal3_Bd_Bv"                         
# [37] "LIB210527RC_AB_1A"                                   
# [38] "LIB210527RC_AB_1B"                                   
# [39] "VC_snRNAseq_12_Animal2_Central_Nucleus"              
# [40] "VC_snRNAseq_12_Animal3_Central_Nucleus"              
# [41] "VC_snRNAseq_12_Animal4_Central_Nucleus"              
# [42] "VC_snRNAseq_12_Animal5_Accessory_Basal__AP2_AP3_"    
# [43] "VC_snRNAseq_12_Animal5_Basal__AP1_"                  
# [44] "VC_snRNAseq_12_Animal5_Basal__AP2_"                  
# [45] "VC_snRNAseq_12_Animal5_Lateral__AP1_"                
# [46] "VC_snRNAseq_12_Animal5_Lateral__AP2_"                
# [47] "VC_snRNAseq_13_Animal5_Central_Nucleus__AP3_" 

unique(metadata$Sample)
# [1] "3c-AMYBLA"                                           
# [2] "4c-AMYBLA"                                           
# [3] "5c-AMYBLA"                                           
# [4] "34ac_scp"                                            
# [5] "35ac_scp"                                            
# [6] "VC_snRNAseq-7_LateralVentralAP1_-1"                  
# [7] "VC_snRNAseq-7_LateralVentralAP2_-2"                  
# [8] "VC_snRNAseq-7_BasalDorsalAP1_-3A"                    
# [9] "VC_snRNAseq-7_BasalDorsalAP1_-3B"                    
# [10] "VC_snRNAseq-7_BasalVentralAP1_-4"                    
# [11] "VC_snRNAseq-7_BasalVentralAP2_-5"                    
# [12] "VC_snRNAseq-7_AccBasalAP1AP2_-7"                     
# [13] "VC_snRNAseq-7_LateralDorsalAP1AP2_-8"                
# [14] "VC_snRNAseq_9_Animal4_LV2"                           
# [15] "VC_snRNAseq_9_Animal4_LV1"                           
# [16] "VC_snRNAseq_9_Animal4_LD"                            
# [17] "VC_snRNAseq_9_Animal4_L-Comb"                        
# [18] "VC_snRNAseq_9_Animal4_BV"                            
# [19] "VC_snRNAseq_9_Animal4_BD"                            
# [20] "VC_snRNAseq_9_Animal4_B-Comb"                        
# [21] "VC_snRNAseq_9_Animal4_AB"                            
# [22] "VC_snRNAseq_8_Animal3_LV"                            
# [23] "VC_snRNAseq_8_Animal3_LV3"                           
# [24] "VC_snRNAseq_8_Animal3_LD"                            
# [25] "VC_snRNAseq_8_Animal3_BV2"                           
# [26] "VC_snRNAseq_8_Animal3_BV1"                           
# [27] "VC_snRNAseq_8_Animal3_BD"                            
# [28] "VC_snRNAseq_8_Animal3_AB__"                          
# [29] "VC_snRNAseq_8_Animal3_Bd_Bv"                         
# [30] "LIB210527RC_AB_1A"                                   
# [31] "LIB210527RC_AB_1B"                                   
# [32] "VC_snRNAseq_12_Animal2_Central_Nucleus"              
# [33] "VC_snRNAseq_12_Animal3_Central_Nucleus"              
# [34] "VC_snRNAseq_12_Animal4_Central_Nucleus"              
# [35] "VC_snRNAseq_12_Animal5_Accessory_Basal__AP2_AP3_"    
# [36] "VC_snRNAseq_12_Animal5_Basal__AP1_"                  
# [37] "VC_snRNAseq_12_Animal5_Basal__AP2_"                  
# [38] "VC_snRNAseq_12_Animal5_Lateral__AP1_"                
# [39] "VC_snRNAseq_12_Animal5_Lateral__AP2_"                
# [40] "VC_snRNAseq_13_Animal5_Central_Nucleus__AP3_"        
# [41] "SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1"         
# [42] "SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1"          
# [43] "SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1"          
# [44] "SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1"          
# [45] "SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1"
# [46] "SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1"          
# [47] "SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1"


# from the metadata, add the "Subregion" and "DV_axis" column to the "sce" sce. Match based on Sample

library("data.table")
sample_info <- metadata


# Create the key column for sce
sce$key <- paste0(sce$Barcode, "_", sce$Sample)

# Convert your data to data.table objects
sce_colData <- as.data.table(colData(sce))
sample_info_dt <- as.data.table(sample_info)

# Identify the common columns to merge on
common_cols <- intersect(names(sce_colData), names(sample_info_dt))

# Perform the merge using the [.data.table method
new_col <- sample_info_dt[sce_colData, on = common_cols]

new_col <- new_col[match(sce$key, sce_colData$key), ]

# Check that the keys match
stopifnot(identical(sce$key, new_col$key))

# Convert new_col to a DataFrame as required by SingleCellExperiment
new_col_df <- DataFrame(new_col)

# Update the colData
colData(sce) <- new_col_df
colnames(sce) <- sce$Barcode

unique(colnames(colData(sce)))
# [1] "Sample_num"            "Sample"                "Species"              
# [4] "Subject"               "Sex"                   "Region"               
# [7] "Subregion"             "DV_axis"               "PI.NeuN"              
# [10] "batch"                 "Barcode"               "sum"                  
# [13] "detected"              "subsets_Mito_sum"      "subsets_Mito_detected"
# [16] "subsets_Mito_percent"  "total"                 "high_mito"            
# [19] "low_lib"               "low_genes"             "discard_auto"         
# [22] "doubletScore"          "sizeFactor"            "species"              
# [25] "key"   

unique(sce$Subregion)

# save sce


# ========= Checking for clusters of doublets ========

png(here(plot_dir, "UMAP_doublet_scores.png"), width=10, height=10, units="in", res=300)
plotReducedDim(sce, dimred = "UMAP", colour_by = "doubletScore")
dev.off()

png(here(plot_dir, "UMAP_species.png"), width=10, height=10, units="in", res=300)
plotReducedDim(sce, dimred = "UMAP", colour_by = "Species")
dev.off()

png(here(plot_dir, "UMAP_ident.png"), width=10, height=10, units="in", res=300)
plotReducedDim(sce, dimred = "UMAP", colour_by = "ident", text_by="ident", point_size=0.1)
dev.off()

# NOTES: 
# It looks like clusters 29, 27, and 19 are mostly doublets. Will remove. 


# ========== plotting UMAPs of broad marker genes ========

# plot some broad marker genes
png(here(plot_dir, "UMAP_broad_markers.png"), width=25, height=25, units="in", res=300)
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

png(here(plot_dir, "Heatmap_fine_markers.png"), width=10, height=10, units="in", res=300)
plotGroupedHeatmap(sce, 
                   group = "ident",
                   features = features,
                   #symmetric=TRUE,
                   center=TRUE,
                   block="Species"
                   )
dev.off()

# dot plots
png(here(plot_dir, "Dotplot_broad_markers.png"), width=10, height=10, units="in", res=300)
plotDots(sce, 
         group = "ident",
         features = features
        )
dev.off()


# ======== Annotation broad cell types ========

excit <- c(0,1,2,8,12,16,20,22,25)
inhib <- c(9,10,11,13,15,17,18,21,24,26,28)
non_neuronal <- c(3,4,5,6,7,14,23,30,31,32,33)
drop <- c(19,27,29)

sce$ident <- as.factor(sce$ident)

# create new broad_celltype column based on above clusters
colData(sce)$broad_celltype <- "NA"
colData(sce)$broad_celltype[which(colData(sce)$ident %in% excit)] <- "Excitatory"
colData(sce)$broad_celltype[which(colData(sce)$ident %in% inhib)] <- "Inhibitory"
colData(sce)$broad_celltype[which(colData(sce)$ident %in% non_neuronal)] <- "Non-neuronal"
colData(sce)$broad_celltype[which(colData(sce)$ident %in% drop)] <- "Drop"

# drop Drop cells
sce <- sce[,sce$broad_celltype != "Drop"]

# plot UMAP with broad cell types
png(here(plot_dir, "UMAP_broad_annotations.png"), width=10, height=10, units="in", res=300)
plotReducedDim(sce, dimred = "UMAP", colour_by = "broad_celltype", point_size=0.1)
dev.off()

# save broad annotations
saveRDS(sce, here("processed-data", "07_annotation", "sce_broad_annotations.rds"))
