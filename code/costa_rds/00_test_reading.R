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

# directories
raw_dir = here("raw-data")
samples_dir = here(raw_dir, "samples")
sampleinfo_dir = here(raw_dir, "sampleinfo")
processed_dir = here("processed-data")


# data paths
barcode.path <- paste0(data_dir, "/Bbn2_Basal1/barcodes.tsv.gz")
features.path <- paste0(data_dir, "/Bbn2_Basal1/features.tsv.gz")
matrix.path <- paste0(data_dir, "/Bbn2_Basal1/matrix.mtx.gz")

mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

mat
# 34702 x 6855 sparse Matrix of class "dgTMatrix"

# Create SCE
sce <- SingleCellExperiment(assays = (list(counts = as(mat,"matrix"))))

# ========= quality control ===========

location <- rowRanges(sce)
is.mito <- any(seqnames(sce) == "MT")

library(scuttle)
df <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
summary(df$subsets_Mito_sum)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 503    3124    5039    5776    7696   38904 

sce <- addPerCellQCMetrics(sce, subsets=list(Mito=is.mito))
colnames(colData(sce))
# [1] "sum"                   "detected"              "subsets_Mito_sum"      "subsets_Mito_detected"
# [5] "subsets_Mito_percent"  "total"

reasons <- perCellQCFilters(df, 
                            sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))
colSums(as.matrix(reasons))


# ====== Normalization ======

clust.sce <- quickCluster(sce) 
sce <- computeSumFactors(sce, cluster=clust.sce, min.mean=0.1)
sce <- logNormCounts(sce)



# Dimensionality reduction
library(scran)
top.genes <- getTopHVGs(sce, n=2000)

set.seed(100) # See below.
sce <- fixedPCA(sce, subset.row=top.genes) 
reducedDimNames(sce)


plotReducedDim(sce, dimred="PCA", colour_by="CORT")
