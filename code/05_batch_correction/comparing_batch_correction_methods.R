# srun -n 1 --mem=64G --cpus-per-task=8 --pty bash -i
library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(batchelor)
library(harmony)
library(Seurat)
library(SeuratData) 
#library(SeuratWrappers)

## save directories
plot_dir = here("plots", "05_batch_correction", "method_comparisons")
processed_dir = here("processed-data","05_batch_correction")


# load sce
load(file=here("processed-data","04_normalization", "sce.human_normalized.rds"))
sce.human
# class: SingleCellExperiment 
# dim: 36601 21268 
# metadata(1): Samples
# assays(1): counts
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(21268): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 5_TTTGTTGTCGGACTGC-1 5_TTTGTTGTCGTTGTTT-1
# colData names(27): Sample Barcode ... discard_auto doubletScore
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

load(file=here("processed-data", "04_normalization", "sce.macaque_normalized.rds"))
sce.macaque
# class: SingleCellExperiment 
# dim: 21369 109142 
# metadata(1): Samples
# assays(1): counts
# rownames(21369): ENSMMUG00000023296 ZNF692 ... ND6 CYTB
# rowData names(7): source type ... gene_biotype Symbol.uniq
# colnames(109142): 1_AAACGAAAGCGCCGTT-1 1_AAACGAACAGCCGTTG-1 ...
# 35_TTTGGTTGTGGCTAGA-1 35_TTTGTTGCACATGGTT-1
# colData names(13): Sample Barcode ... discard_auto doubletScore
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

load(file=here("processed-data","04_normalization", "sce.baboon_normalized.rds"))
sce.baboon
# class: SingleCellExperiment 
# dim: 20842 57564 
# metadata(1): Samples
# assays(1): counts
# rownames(20842): PITHD1 ACTL8 ... ENSPANG00000037790 ENSPANG00000045064
# rowData names(7): source type ... gene_biotype Symbol.uniq
# colnames(57564): 1_AAACCCAAGACTGTTC-1 1_AAACCCAAGATCCAAA-1 ...
# 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# colData names(21): Sample Barcode ... discard_auto doubletScore
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):


# ===== Feature seleection =====

dec1 <- modelGeneVar(sce.human)
dec2 <- modelGeneVar(sce.macaque)
dec3 <- modelGeneVar(sce.baboon)

combined.dec <- combineVar(dec1, dec2, dec3)
chosen.hvgs <- getTopHVGs(combined.dec, n=5000)

# As a sanity check, let's confirm that there are indeed batch effects in the data

combined <- correctExperiments(Human=sce.human, Macaque=sce.macaque, Baboon=sce.baboon, PARAM=NoCorrectParam())
combined
# class: SingleCellExperiment 
# dim: 14391 187974 
# metadata(0):
#   assays(3): merged counts logcounts
# rownames(14391): SAMD11 NOC2L ... EIF1AY RPS4Y2
# rowData names(4): gene_name gene_type Symbol.uniq gene_biotype
# colnames(187974): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ... 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# colData names(15): batch Sample ... doubletScore sizeFactor
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):



# Initialize an empty list to store subsets
subsets <- list()

# Loop through each batch
for (batch_id in unique(combined$batch)) {
    # Filter cells from the current batch
    batch_cells <- combined[, combined$batch == batch_id]
    
    # Randomly sample 6000 nuclei from this batch
    set.seed(123) # for reproducibility
    sampled_cells <- batch_cells[, sample(ncol(batch_cells), 6000)]
    
    # Add the sampled cells to the list
    subsets[[batch_id]] <- sampled_cells
}

# Combine the subsets from each batch
combined.subset <- do.call(cbind, subsets)
combined.subset
# class: SingleCellExperiment 
# dim: 14391 18000 
# metadata(0):
#   assays(3): merged counts logcounts
# rownames(14391): SAMD11 NOC2L ... EIF1AY RPS4Y2
# rowData names(4): gene_name gene_type Symbol.uniq gene_biotype
# colnames(18000): 5_GACCCTTCATGCAGGA-1 5_GACTGATAGGGACACT-1 ... 4_CTGAATGTCTTCCTAA-1 7_CGTCCATTCCTCGCAT-1
# colData names(15): batch Sample ... doubletScore sizeFactor
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):

rm(combined)


# ============== PCA and TSNE before correction ==================

combined.subset$Species <- combined.subset$batch


set.seed(111)
combined.subset <- runPCA(combined.subset, subset_row=chosen.hvgs,
                          BPPARAM = BiocParallel::MulticoreParam(workers=8))

combined.subset <- runTSNE(combined.subset, dimred="PCA",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))

pdf(here(plot_dir, "TSNE_batch_effects.pdf"))
plotTSNE(combined.subset, colour_by="batch")
dev.off()


# ============== MNN on Species > Samples ============


# Run batch correction with MNN across Samples
mnn_species_pca <- batchelor::reducedMNN(reducedDim(combined.subset, "PCA"),
                                         batch=as.factor(combined.subset$Species))

reducedDim(combined.subset,"mnn_species_pca") <- mnn_species_pca$corrected


# Run batch correction with MNN across Species 
mnn_species_sample_pca <- batchelor::reducedMNN(reducedDim(combined.subset, "mnn_species_pca"),
                                                batch=as.factor(combined.subset$Sample))

reducedDim(combined.subset,"mnn_species_sample_pca") <- mnn_species_sample_pca$corrected


# compute TSNE
combined.subset <- runTSNE(combined.subset, dimred="mnn_species_sample_pca",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))


# plot
pdf(here(plot_dir, "TSNE_mnn_species_sample.pdf"))
plotReducedDim(combined.subset, dim="TSNE", colour_by="Species")
dev.off()



# ================ Hrmony on Species > Samples ==============

combined.subset <- RunHarmony(combined.subset, c("Species", "Sample"))
combined.subset

# compute TSNE
combined.subset <- runTSNE(combined.subset, dimred="HARMONY",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))


# plot
pdf(here(plot_dir, "TSNE_harmony_species_samples.pdf"))
plotReducedDim(combined.subset, dim="TSNE", colour_by="Species")
dev.off()


# ================ Hrmony on Samples > Species ==============

combined.subset <- RunHarmony(combined.subset, c("Species"))
combined.subset

# compute TSNE
combined.subset <- runTSNE(combined.subset, dimred="HARMONY",
                           BPPARAM = BiocParallel::MulticoreParam(workers=8))


# plot
pdf(here(plot_dir, "TSNE_harmony_species.pdf"))
plotReducedDim(combined.subset, dim="TSNE", colour_by="Species")
dev.off()


# ================ Seurat ==============
# remove duplicate columns 
combined.subset <- combined.subset[, !duplicated(colnames(combined.subset))]

combined.seurat <- as.Seurat(combined.subset, counts="counts", data="logcounts")
combined.seurat
# An object of class Seurat 
# 14391 features across 18000 samples within 1 assay 
# Active assay: originalexp (14391 features, 0 variable features)
# 2 layers present: counts, data
# 6 dimensional reductions calculated: PCA, TSNE, mnn_species_pca, mnn_samples_pca, mnn_species_sample_pca, mnn_samples_species_pca

DimPlot(combined.seurat, reduction="TSNE", group.by="Species")

combined.seurat <- NormalizeData(combined.seurat)
combined.seurat <- FindVariableFeatures(combined.seurat)
combined.seurat <- ScaleData(combined.seurat)

#combined.seurst <- SCTransform(combined.seurat, verbose = FALSE)
combined.seurat <- RunPCA(combined.seurat)

combined.seurat <- FindNeighbors(combined.seurat, dims=1:30, reduction="pca")
combined.seurat <- FindClusters(combined.seurat, resolution=2, cluster.name = "unintegrated_clusters")
combined.seurat <- RunUMAP(combined.seurat, dims=1:30, reduction="pca", reduction.name = "umap.unintegrated")


DimPlot(combined.seurat, reduction = "umap.unintegrated", group.by = c("Species"))


# split the dataset into a list of two seurat objects (stim and CTRL)
seurat.list <- SplitObject(combined.seurat, split.by = "Species")

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

seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features, reduction = "rpca")

# this command creates an 'integrated' data assay
seurat.int <- IntegrateData(anchorset = seurat.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(seurat.int) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.int <- ScaleData(seurat.int, verbose = FALSE)
seurat.int <- RunPCA(seurat.int, npcs = 30, verbose = FALSE)
seurat.int <- RunTSNE(seurat.int, reduction = "pca", dims = 1:30)
seurat.int <- FindNeighbors(seurat.int, reduction = "pca", dims = 1:30)
seurat.int <- FindClusters(seurat.int, resolution = 0.5)

# Visualization
pdf(here(plot_dir, "Seurat_UMAP_unintegrated.pdf"))
DimPlot(seurat.int, reduction = "tsne", group.by = "Species")
dev.off()


