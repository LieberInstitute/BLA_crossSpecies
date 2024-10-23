library("spatialLIBD")
library("patchwork")
library("scran")
library("here")
library("scater")

# Save directories
processed_dir <- here("processed-data", "09_psychiatric_GSEA")

# load sce
sce.excit <- readRDS(here("processed-data", "07_annotation_and_characterization", "sce_excit_final_subclusters_annotated.rds"))
sce.inhib <- readRDS(here("processed-data", "07_annotation_and_characterization", "sce_inhib_final_subclusters_annotated.rds"))
sce.other <- readRDS(here("processed-data", "07_annotation_and_characterization", "sce_other_final_celltypes.rds"))

# combine SCE objects
sce.inhib$subcluster_idents <- NULL
sce.inhib$integrated_snn_res.0.4 <- NULL
sce.excit$integrated_snn_res.0.4 <- NULL

sce.neurons <- cbind(sce.excit, sce.inhib)
sce <- cbind(sce.neurons, sce.other)

# remove objects to save memory
rm(sce.excit, sce.inhib, sce.other, sce.neurons)
gc()


# subset to just humans
sce.human <- sce[, sce$species == "human"]
sce.human
# class: SingleCellExperiment 
# dim: 13842 20625 
# metadata(0):
# assays(2): counts logcounts
# rownames(13842): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(0):
# colnames(20625): AAACCCACAGGTCCCA-1 AAACCCATCCATTGTT-1 ...
#   TTTGGTTCAGCACACC-1 TTTGTTGTCGGACTGC-1
# colData names(34): orig.ident nCount_originalexp ... broad_celltype
#   fine_celltype
# reducedDimNames(2): PCA UMAP
# mainExpName: originalexp
# altExpNames(0):


# get percent neurons to non-neurons in "broad_celltype"
table(sce.human$broad_celltype)

# ======= Load original human data =======

sce.human.og <- readRDS(here("processed-data", "04_normalization", "sce.human_normalized.rds"))
sce.human.og
# class: SingleCellExperiment 
# dim: 21844 21212 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(21844): AL627309.1 AL627309.5 ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(21212): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
#   5_TTTGTTGTCGGACTGC-1 5_TTTGTTGTCGTTGTTT-1
# colData names(31): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):




# ======== Copying over meta data =========

# we can see that ~600 cells have been dropped from the original data. So we will need to only retain the matched cells
# and then copy the metadata over. 

# get the cells that are in both datasets
sce.human.new <- sce.human.og[,sce.human.og$key %in% sce.human$key]
sce.human.new
# class: SingleCellExperiment 
# dim: 21844 20625 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(21844): AL627309.1 AL627309.5 ... AC007325.4 AC007325.2
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(20625): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
#   5_TTTGTTGTCGGACTGC-1 5_TTTGTTGTCGTTGTTT-1
# colData names(31): Sample Barcode ... doubletScore sizeFactor
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):

# get colData
og_meta <- data.frame(colData(sce.human.new))
final_meta <- data.frame(colData(sce.human))

# Reorder final_metadata to match the order of keys in og_metadata
final_meta_ordered <- final_meta[match(og_meta$key, final_meta$key), ]
identical(og_meta$key, final_meta_ordered$key)  
# TRUE

# copy over metadata
colData(sce.human.new) <- as(final_meta_ordered, "DataFrame")
sce <- sce.human.new

# save SCE
saveRDS(sce, here("processed-data", "09_psychiatric_GSEA", "sce.human_all_genes.rds"))

#load SCE
sce <- readRDS(here("processed-data", "09_psychiatric_GSEA", "sce.human_all_genes.rds"))


# ======= Getting human cell type marker genes for deconvolution ========
library(DeconvoBuddies)
library(dplyr)

# drop mito genes
sce <- sce[!grepl("^MT-", rownames(sce)),]

# filter cell types with < 70 cells
celltype_counts <- table(colData(sce)$fine_celltype)
celltypes_to_drop <- names(celltype_counts[celltype_counts < 70])
sce <- sce[, !(colData(sce)$fine_celltype %in% celltypes_to_drop)]
sce$fine_celltype <- factor(sce$fine_celltype)

# make sure the data are properly normalized
clust <- quickCluster(sce) 
sce <- computeSumFactors(sce, cluster=clust, min.mean=0.1)
sce <- logNormCounts(sce)



# 1 vs All
oneVall <- findMarkers_1vAll(
    sce,
    assay_name = "logcounts",
    cellType_col = "fine_celltype",
    add_symbol = FALSE,
    mod = "~Sample",
    verbose = TRUE
)

top_markers <- oneVall %>%
    group_by(cellType.target) %>%
    top_n(n = 100, wt = logFC)

write.csv(top_markers, here(processed_dir, "top_100_human_markers_fine_celltype.csv"), row.names = FALSE)

# mean ratio markers
ratio_markers <- get_mean_ratio(
    sce,
    cellType_col = "fine_celltype",
    assay_name = "logcounts"
)

top_ratio_markers <- ratio_markers %>%
    group_by(cellType.target) %>%
    top_n(n = 100, wt = MeanRatio)

write.csv(top_ratio_markers, here(processed_dir, "top_100_human_ratio_markers_fine_celltype.csv"), row.names = FALSE)


