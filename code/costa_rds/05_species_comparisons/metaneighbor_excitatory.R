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



# ====== Subset inhibitory celltypes ======

sce.totty <- sce.totty[, sce.totty$celltype == "Excit"]
sce.yu <- sce.yu[, sce.yu$celltype == "ExN"]


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


# ========== Identify top hits ==========

top_hits = topHits(cell_NV = celltype_NV,
                   dat = sce.bla,
                   study_id = sce.bla$study_id,
                   cell_type = sce.bla$cell_type,
                   threshold = 0.9)
top_hits
# Study_ID|Celltype_1 Study_ID|Celltype_2 Mean_AUROC         Match_type
# 1         Yu|LAMP5 ABO        Totty|Exc_05       0.98 Reciprocal_top_hit
# 2        Yu|SATB2 IL15        Totty|Exc_10       0.98 Reciprocal_top_hit
# 3        Yu|SATB2 IL15        Totty|Exc_07       0.98          Above_0.9
# 4         Totty|Exc_10    Yu|SATB2 ST8SIA2       0.97          Above_0.9
# 5        Yu|SATB2 IL15        Totty|Exc_08       0.97          Above_0.9
# 6         Totty|Exc_08           Yu|STRIP2       0.97          Above_0.9
# 7         Totty|Exc_07     Yu|SATB2 CALCRL       0.96          Above_0.9
# 8         Totty|Exc_08           Yu|TFAP2C       0.95          Above_0.9
# 9      Yu|HGF C11orf87        Totty|Exc_04       0.94 Reciprocal_top_hit
# 10    Yu|LAMP5 COL25A1        Totty|Exc_03       0.94 Reciprocal_top_hit
# 11       Yu|SATB2 IL15        Totty|Exc_11       0.94          Above_0.9
# 12        Totty|Exc_04        Yu|HGF NPSR1       0.92          Above_0.9
# 13      Yu|RXFP2 RSPO2        Totty|Exc_02       0.92 Reciprocal_top_hit
# 14        Totty|Exc_03       Yu|LAMP5 BDNF       0.91          Above_0.9
# 15        Totty|Exc_08       Yu|SOX11 EBF2       0.91          Above_0.9
# 16       Yu|VGLL3 MEPE        Totty|Exc_12       0.85 Reciprocal_top_hit