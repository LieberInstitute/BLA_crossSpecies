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
plot_dir = here("plots", "05_batch_correction", "seurat_v4")
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

# ========= Conver to Seurat v4 ==========

combined$species <- combined$batch

# rename $species. 1 = human, 2 = baboon, 3 = macaque
combined$species[combined$species == 1] <- "human"
combined$species[combined$species == 2] <- "baboon"
combined$species[combined$species == 3] <- "macaque"
unique(combined$species)
# [1] "human"   "baboon"  "macaque"

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
save(seurat.anchors, file = here(processed_dir, "seurat.anchors.rda"))

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
save(seurat.int, file = here(processed_dir, "seurat_integrated.rda"))


# ======= high quality plotting =======
load(seurat.int, file = here(processed_dir, "seurat_integrated.rda"))

pdf(here(plot_dir, "UMAP_species_integrated_by_samples.pdf"))
DimPlot(seurat.int, reduction = "umap", group.by = "species", raster=FALSE)
dev.off()

pdf(here(plot_dir, "UMAP_samples_integrated_by_samples.pdf"))
DimPlot(seurat.int, reduction = "umap", group.by = "Sample", , raster=FALSE) + NoLegend()
dev.off()

