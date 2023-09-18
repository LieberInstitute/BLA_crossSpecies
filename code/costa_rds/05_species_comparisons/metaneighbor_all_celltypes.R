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

# Bab.200 <-readRDS(here(processed_dir, "00_costa_rds", "AllBaboon_Amygdala200_04.rds"))
# Bab.200
# sce.bab <- as.SingleCellExperiment(Bab.200)

human.yu <- readRDS(here(processed_dir, "00_costa_rds", "GSE195445_Human_obj.rda"))
human.yu
sce.yu <- as.SingleCellExperiment(human.yu)

# remove "Human" form the beginning of $ident
levels(sce.yu$ident) <- gsub("Human_", "", levels(sce.yu$ident))

load(here(processed_dir, "00_costa_rds", "sce_annotated.rda"))
sce.totty <- sce


# ====== MetaNeighbor ======

common_genes <- intersect(rownames(sce.totty), rownames(sce.yu))
sce.totty <- sce.totty[common_genes,]
sce.yu <- sce.yu[common_genes, ]

new_colData = data.frame(
  study_id = rep(c('Yu', 'Totty'), c(ncol(sce.yu), ncol(sce.totty))),
  cell_type = c(as.character(colData(sce.yu)$ident), colData(sce.totty)$annotation)
)
sce.bla <- SingleCellExperiment(
  Matrix(cbind(assay(sce.yu, 1), assay(sce.totty, 1)), sparse = TRUE),
  colData = new_colData
)
dim(sce.bla)
rm(sce.totty); rm(sce.yu)

# ====== run it ======

var_genes = variableGenes(dat = sce.bla, exp_labels = sce.bla$study_id)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = sce.bla,
                             study_id = sce.bla$study_id,
                             cell_type = sce.bla$cell_type,
                             fast_version = TRUE)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
gplots::heatmap.2(celltype_NV,
                  margins=c(8,8),
                  keysize=1,
                  key.xlab="AUROC",
                  key.title="NULL",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.7,
                  cexCol = 0.7)
