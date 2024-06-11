# srun -n 1 --mem=400G --cpus-per-task=8 --pty bash -i
#remotes::install_version("Matrix", version = "1.6-1")
library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(batchelor)
library(harmony)
library(Seurat)
#library(future)


## save directories
plot_dir = here("plots", "05_batch_correction", "seurat_v4", "without_CeA")
processed_dir = here("processed-data","05_batch_correction")

# Enable parallelization
#plan("multicore", workers = 10)

# load sce
sce.human <- readRDS(file=here("processed-data","04_normalization", "sce.human_normalized.rds"))
sce.human
# class: SingleCellExperiment 
# dim: 13874 21268 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(21268): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 5_TTTGTTGTCGGACTGC-1 5_TTTGTTGTCGTTGTTT-1
# colData names(28): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

sce.macaque <- readRDS(file=here("processed-data", "04_normalization", "sce.macaque_normalized.rds"))
sce.macaque
# class: SingleCellExperiment 
# dim: 13874 109142 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(7): source type ... gene_biotype Symbol.uniq
# colnames(109142): 1_AAACGAAAGCGCCGTT-1 1_AAACGAACAGCCGTTG-1 ...
# 35_TTTGGTTGTGGCTAGA-1 35_TTTGTTGCACATGGTT-1
# colData names(14): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

sce.baboon <- readRDS(file=here("processed-data","04_normalization", "sce.baboon_normalized.rds"))
sce.baboon
# class: SingleCellExperiment 
# dim: 13874 47039 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(5): ID Symbol Type gene_id_ncbi Symbol.uniq
# colnames(47039): 1_AAACCCAAGACTGTTC-1 1_AAACCCAAGATCCAAA-1 ...
# 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# colData names(22): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):



 # ====== Subset to only 1:1 orthologs ======

 method <- "gprofiler"

 baboon_orthos <- orthogene::convert_orthologs(gene_df = counts(sce.baboon),
                                         gene_input = "rownames", 
                                         gene_output = "rownames", 
                                         input_species = "macaque",
                                         output_species = "human",
                                         non121_strategy = "drop_both_species",
                                         method = method) 
 # =========== REPORT SUMMARY ===========
 # Total genes dropped after convert_orthologs :
 #     5,487 / 21,091 (26%)
 # Total genes remaining after convert_orthologs :
 #     15,604 / 21,091 (74%)
 macaque_orthos <- orthogene::convert_orthologs(gene_df = counts(sce.macaque),
                                               gene_input = "rownames", 
                                               gene_output = "rownames", 
                                               input_species = "baboon",
                                               output_species = "human",
                                               non121_strategy = "drop_both_species",
                                               method = method) 
 # =========== REPORT SUMMARY ==========
 # Total genes dropped after convert_orthologs :
 #     6,221 / 21,369 (29%)#
# Total genes remaining after convert_orthologs :
#     15,148 / 21,369 (71%)

 nhp_orthologs <- intersect(rownames(baboon_orthos), rownames(macaque_orthos))
 length(nhp_orthologs)
 # [1] 14025

 # get only orthologs that are in the human dataset
 valid.human.orthologs <- nhp_orthologs[nhp_orthologs %in% rownames(sce.human)]
 length(valid.human.orthologs)
 # [1] 13948

 # get only orthologs that are in the macaque dataset
 valid.macaque.orthologs <- nhp_orthologs[nhp_orthologs %in% rownames(sce.macaque)]
 length(valid.macaque.orthologs)
 # [1] 13952

 # get only orthologs that are in the baboon dataset
 valid.baboon.orthologs <- nhp_orthologs[nhp_orthologs %in% rownames(sce.baboon)]
 length(valid.baboon.orthologs)
 # [1] 13937

 # get only orthologs common to all datasets
 human_v_mac.orthologs <- intersect(valid.human.orthologs, valid.macaque.orthologs)
 final.121.orthologs <- intersect(human_v_mac.orthologs, valid.baboon.orthologs)
 length(final.121.orthologs)
 # [1] 13874

 # subset to only 1:1 orthologs
 sce.human.ortho <- sce.human[final.121.orthologs ,]
 sce.human.ortho
 # class: SingleCellExperiment 
 # dim: 13874 21268 
 # metadata(1): Samples
 # assays(1): counts
 # rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
 # rowData names(7): source type ... gene_type Symbol.uniq
 # colnames(21268): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
 # 5_TTTGTTGTCGGACTGC-1 5_TTTGTTGTCGTTGTTT-1
 # colData names(27): Sample Barcode ... discard_auto doubletScore
 # reducedDimNames(0):
 #     mainExpName: NULL
 # altExpNames(0):

 sce.macaque.ortho <- sce.macaque[final.121.orthologs,]
 sce.macaque.ortho
 # dim: 13874 109142 
 # metadata(1): Samples
 # assays(1): counts
 # rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
 # rowData names(7): source type ... gene_biotype Symbol.uniq
 # colnames(109142): 1_AAACGAAAGCGCCGTT-1 1_AAACGAACAGCCGTTG-1 ...
 # 35_TTTGGTTGTGGCTAGA-1 35_TTTGTTGCACATGGTT-1
 # colData names(13): Sample Barcode ... discard_auto doubletScore
 # reducedDimNames(0):
 #     mainExpName: NULL
 # altExpNames(0):

 sce.baboon.ortho <- sce.baboon[final.121.orthologs,]
 sce.baboon.ortho
 # class: SingleCellExperiments
 # dim: 13874 47039 
 # metadata(1): Samples
 # assays(1): counts
 # rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
 # rowData names(5): ID Symbol Type gene_id_ncbi Symbol.uniq
 # colnames(47039): 1_AAACCCAAGACTGTTC-1 1_AAACCCAAGATCCAAA-1 ...
 # 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
 # colData names(21): Sample Barcode ... discard_auto doubletScore
 # reducedDimNames(0):
 #     mainExpName: NULL
 # altExpNames(0):

 # rename variables
sce.human <- sce.human.ortho
sce.macaque <- sce.macaque.ortho
sce.baboon <- sce.baboon.ortho

dim(sce.human)
# [1] 13874 21268

dim(sce.macaque)
# [1] 13874 109142

dim(sce.baboon)
# [1] 13874 47039


# ========== Combine datasets without correction using Batchelor ==========

combined <- correctExperiments(sce.human, sce.baboon, sce.macaque, 
                               PARAM=NoCorrectParam())


# ============ Batch correction with Seurat v4 (RPCA) ============
# remove duplicate columns 
combined<- combined[, !duplicated(colnames(combined))]
combined
# dim: 14391 187856 
# metadata(0):
#     assays(4): merged counts logcounts binomial_deviance_residuals
# rownames(14391): SAMD11 NOC2L ... EIF1AY RPS4Y2
# rowData names(4): gene_name gene_type Symbol.uniq gene_biotype
# colnames(187856): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# colData names(16): batch Sample ... sizeFactor Species
# reducedDimNames(4): PCA TSNE GLM-PCA_approx UMAP
# mainExpName: NULL
# altExpNames(0):

colnames(colData(combined))
# [1] "batch"                 "Sample"                "Barcode"              
# [4] "sum"                   "detected"              "subsets_Mito_sum"     
# [7] "subsets_Mito_detected" "subsets_Mito_percent"  "total"                
# [10] "high_mito"             "low_lib"               "low_genes"            
# [13] "discard_auto"          "doubletScore"          "sizeFactor" 


# save RAM by removing the original datasets
rm(sce.human)
rm(sce.baboon)
rm(sce.macaque)

combined$species <- combined$batch

# rename $species. 1 = human, 2 = baboon, 3 = macaque
combined$species[combined$species == 1] <- "human"
combined$species[combined$species == 2] <- "baboon"
combined$species[combined$species == 3] <- "macaque"
unique(combined$species)
# [1] "human"   "baboon"  "macaque"



# Read in the CSV file as a data frame
tmp <- read.delim(here("raw-data","sampleinfo",
                       "master_sampleinfo_2023-10-09.csv"),
                  header = T,sep=',')

# View the subsetted data
head(tmp)


sce <- combined
rm(combined)


#        Sample Species Subject Sex   Region Subregion DV_axis  PI.NeuN
#   1 3c-AMYBLA   Human  Br2327     Amygdala       BLA         PI+NeuN+
#   2 4c-AMYBLA   Human  Br8692     Amygdala       BLA         PI+NeuN+
#   3 5c-AMYBLA   Human  Br9021     Amygdala       BLA         PI+NeuN+
#   4  34ac_scp   Human  Br8331     Amygdala       BLA         PI+NeuN+
#   5  35ac_scp   Human  Br5273     Amygdala       BLA         PI+NeuN+
#   6   Sample1 Macaque    Mac2     Amygdala   Lateral Ventral         

# get only amygdala region
metadata <- tmp[tmp$Region == "Amygdala",]

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

combined <- sce
rm(sce)
save(combined, file = here(processed_dir, "sce_combined.rda"))


# drop Central Nucleus samples
combined <- combined[, !grepl("Central Nucleus", colnames(combined))]





# ========= Conver to Seurat v4 ==========

combined.seurat <- as.Seurat(combined, counts="counts", data="logcounts")
combined.seurat
# An object of class Seurat 
# 13874 features across 177346 samples within 1 assay 
# Active assay: originalexp (13874 features, 0 variable features)
# 2 layers present: counts, data

# standard normalization in seurat
combined.seurat <- NormalizeData(combined.seurat)
combined.seurat <- FindVariableFeatures(combined.seurat)
combined.seurat <- ScaleData(combined.seurat)

# Run PCA
combined.seurat <- RunPCA(combined.seurat)

# Run UMAP
combined.seurat <- FindNeighbors(combined.seurat, dims=1:30, reduction="pca")
combined.seurat <- FindClusters(combined.seurat, resolution=2, cluster.name = "unintegrated_clusters")
combined.seurat <- RunUMAP(combined.seurat, dims=1:30, reduction="pca", reduction.name = "umap.unintegrated")

# plot uncorrected data
pdf(here(plot_dir, "UMAP_species_uncorrected.pdf"))
DimPlot(combined.seurat, reduction = "umap.unintegrated", group.by = c("species")) 
dev.off()

pdf(here(plot_dir, "UMAP_samples_uncorrected.pdf"))
DimPlot(combined.seurat, reduction = "umap.unintegrated", group.by = c("Sample")) + NoLegend()
dev.off()



# split the dataset into a list of seurat objects by species
seurat.list <- SplitObject(combined.seurat, split.by = "Sample")

# normalize and identify variable features for each dataset independently
seurat.list<- lapply(X = seurat.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = seurat.list)
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# find integration anchors
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features, reduction = "rpca")

# save integrated seurat object
save(seurat.anchors, file = here(processed_dir, "seurat.anchors_noCeA.rda"))

# integrate the datasets
seurat.int <- IntegrateData(anchorset = seurat.anchors)
DefaultAssay(seurat.int) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.int <- ScaleData(seurat.int, verbose = FALSE)
seurat.int <- RunPCA(seurat.int, npcs = 30, verbose = FALSE)
seurat.int <- RunUMAP(seurat.int, reduction = "pca", dims = 1:30)
seurat.int <- FindNeighbors(seurat.int, reduction = "pca", dims = 1:30)
seurat.int <- FindClusters(seurat.int, resolution = 0.5)

# Visualization
pdf(here(plot_dir, "UMAP_species_integrated_by_samples.pdf"))
DimPlot(seurat.int, reduction = "umap", group.by = "species")
dev.off()

pdf(here(plot_dir, "UMAP_samples_integrated_by_samples.pdf"))
DimPlot(seurat.int, reduction = "umap", group.by = "Sample") + NoLegend()
dev.off()



# save integrated seurat object
save(seurat.int, file = here(processed_dir, "seurat_integrated_noCeA.rda"))


# ======= high quality plotting =======
load(seurat.int, file = here(processed_dir, "seurat_integrated.rda"))

pdf(here(plot_dir, "UMAP_species_integrated_by_samples.pdf"))
DimPlot(seurat.int, reduction = "umap", group.by = "species", raster=FALSE)
dev.off()

pdf(here(plot_dir, "UMAP_samples_integrated_by_samples.pdf"))
DimPlot(seurat.int, reduction = "umap", group.by = "Sample", , raster=FALSE) + NoLegend()
dev.off()

