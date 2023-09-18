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

# directories
raw_dir = here("raw-data")
samples_dir = here(raw_dir, "samples")
sampleinfo_dir = here(raw_dir, "sampleinfo")
processed_dir = here("processed-data")
plot_dir = here("plots", "05_species_comparisons")

Amy.700 <- readRDS(here(processed_dir, "00_costa_rds", "amygdalaCNScaled700.rds"))
Amy.700
sce.costa <- as.SingleCellExperiment(Amy.700)
# An object of class Seurat 
# 21353 features across 120794 samples within 2 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

Bab.200 <-readRDS(here(processed_dir, "00_costa_rds", "AllBaboon_Amygdala200_04.rds"))
Bab.200
sce.bab <- as.SingleCellExperiment(Bab.200)

human.yu <- readRDS(here(processed_dir, "00_costa_rds", "GSE195445_Human_obj.rda"))
human.yu
sce.yu <- as.SingleCellExperiment(human.yu)

# remove "Human" form the beginning of $ident
levels(sce.yu$ident) <- gsub("Human_", "", levels(sce.yu$ident))

load(here(processed_dir, "00_costa_rds", "sce_annotated.rda"))
sce.totty <- sce


# ====== MetaNeighbor ======

common_genes <- intersect(rownames(sce.totty), rownames(sce.yu))

sce.totty <- sce.totty[common_genes, ]
sce.yu <- sce.yu[common_genes, ]


# Combine expression matrices
data <- cbind(assay(sce.totty), assay(sce.yu))

# Create experiment labels
exp.lab <- c(rep("totty", ncol(sce.totty)), rep("yu", ncol(sce.yu)))

# Assuming you have cell type labels in a column named "cell_type" in colData
cell.lab <- c(colData(sce.totty)$celltype, colData(sce.yu)$ident)

# run metaneighbor
AUROC.scores <- MetaNeighbor(dat = data, experiment_labels = exp.lab, celltype_labels = cell.lab, genesets = genesets, bplot=TRUE)

# plot
hist(AUROC.scores, main="Distribution of AUROC Scores", xlab="AUROC Scores", breaks=20, xlim=c(0,1))
abline(v=mean(AUROC.scores), col="red", lty=2, lwd=2)

