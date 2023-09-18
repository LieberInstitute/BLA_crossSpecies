library("dplyr")
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
library("pheatmap")

# load directories
raw_dir = here("raw-data")
samples_dir = here(raw_dir, "samples")
sampleinfo_dir = here(raw_dir, "sampleinfo")

# save directories
processed_dir = here("processed-data","04_anntation")
plot_dir = here("plots", "04_annotation")

# load costa data
mac.700 <- readRDS(here("processed-data", "00_costa_rds", "amygdalaCNScaled700.rds"))
mac.700
# An object of class Seurat 
# 21353 features across 120794 samples within 2 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

bab.200 <- readRDS(here("processed-data", "00_costa_rds", "AllBaboon_Amygdala200_04.rds"))
bab.200
# An object of class Seurat 
# 31544 features across 71903 samples within 2 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap



# ======= convert to SCE ========

sce.mac <- as.SingleCellExperiment(mac.700)
sce.bab <- as.SingleCellExperiment(bab.200)


# ======= UMAP plots =======
features = c("seurat_clusters", "SNAP25", "SLC17A7", "GAD2", "Region", "Sample")

# macaque
pdf(here(plot_dir,"macaque_UMAP_features.pdf"))
for (feature in features) {
  p <- plotUMAP(sce.mac, color_by=feature, text_by="seurat_clusters")
  print(p)
}
dev.off()

# baboon
pdf(here(plot_dir,"baboon_UMAP_features.pdf"))
for (feature in features) {
  p <- plotUMAP(sce.bab , color_by=feature, text_by="seurat_clusters")
  print(p)
}
dev.off()



# ============ Broad Marker Genes ==========

# ==== Dot plots ===== 
broad_markers = c('SNAP25', #neurons
                  'SLC17A7', #excit
                  'SLC17A6', 
                  'GAD1',  #inhib
                  'GAD2', 
                  'AQP4', #astro
                  'MOBP') #WM


# baboon
pdf(width=7, height=7, here(plot_dir,"broadMarkers_expression_baboon.pdf"))
plotExpression(sce.bab, features=broad_markers, x="seurat_clusters", color_by="seurat_clusters", 
               show_median = TRUE, ncol=1)
dev.off()


# macaque
pdf(width=7, height=7, here(plot_dir,"broadMarkers_expression_macaque.pdf"))
plotExpression(sce.mac, features=broad_markers, x="seurat_clusters", color_by="seurat_clusters", 
               show_median = TRUE, ncol=1)
dev.off()


# ===== UMAPs =====

# baboon
pdf(here(plot_dir,"broadMarkers_UMAP_baboon.pdf"))
for (marker in broad_markers) {
  p <- plotUMAP(sce.bab , color_by=marker, text_by="seurat_clusters")
  print(p)
}
dev.off()

# macaque
pdf(here(plot_dir,"broadMarkers_UMAP_macaque.pdf"))
for (marker in broad_markers) {
  p <- plotUMAP(sce.mac , color_by=marker, text_by="seurat_clusters")
  print(p)
}
dev.off()


# ====Non-neuron subtype marker genes ===== 
subtype_markers = c('SNAP25',
                  'AQP4',
                  'MOBP',
                  'CTSS',
                  'CSPG4',
                  'TIE1')
# baboon
pdf(width=7, height=7, here(plot_dir,"subtypeMarkers_non-neurons_expression_baboon.pdf"))
plotExpression(sce.bab, features=subtype_markers, x="seurat_clusters", color_by="seurat_clusters", 
               show_median = TRUE, ncol=1)
dev.off()


# macaque
pdf(width=7, height=7, here(plot_dir,"subtypeMarkers_non-neurons_expression_macaque.pdf"))
plotExpression(sce.mac, features=subtype_markers, x="seurat_clusters", color_by="seurat_clusters", 
               show_median = TRUE, ncol=1)
dev.off()


# ===== UMAPs =====

# baboon
pdf(here(plot_dir,"subtypeMarkers_non-neurons_UMAP_baboon.pdf"))
for (marker in subtype_markers) {
  p <- plotUMAP(sce.bab , color_by=marker, text_by="seurat_clusters")
  print(p)
}
dev.off()

# macaque
pdf(here(plot_dir,"subtypeMarkers_non-neurons_UMAP_macaque.pdf"))
for (marker in subtype_markers) {
  p <- plotUMAP(sce.mac , color_by=marker, text_by="seurat_clusters")
  print(p)
}
dev.off()



# ==== Inhibitory subtype marker genes ===== 
subtype_markers = c('SNAP25',
                    'GAD1',
                    'GAD2',
                    'SST',
                    'VIP',
                    'PVALB',
                    'LAMP5',
                    'FOXP2')
# baboon
pdf(width=7, height=7, here(plot_dir,"subtypeMarkers_inhibitory_expression_baboon.pdf"))
plotExpression(sce.bab, features=subtype_markers, x="seurat_clusters", color_by="seurat_clusters", 
               show_median = TRUE, ncol=1)
dev.off()


# macaque
pdf(width=7, height=7, here(plot_dir,"subtypeMarkers_inhibitory_expression_macaque.pdf"))
plotExpression(sce.mac, features=subtype_markers, x="seurat_clusters", color_by="seurat_clusters", 
               show_median = TRUE, ncol=1)
dev.off()


# ===== UMAPs =====

# baboon
pdf(here(plot_dir,"subtypeMarkers_inhibitory_UMAP_baboon.pdf"))
for (marker in subtype_markers) {
  p <- plotUMAP(sce.bab , color_by=marker, text_by="seurat_clusters")
  print(p)
}
dev.off()

# macaque
pdf(here(plot_dir,"subtypeMarkers_inhibitory_UMAP_macaque.pdf"))
for (marker in subtype_markers) {
  p <- plotUMAP(sce.mac , color_by=marker, text_by="seurat_clusters")
  print(p)
}
dev.off()