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
library("harmony")

## Load sce
load(here("processed-data", "01_build_sce","sce_macaque_raw.rda"))
dim(sce)

# Save directories
plot_dir = here("plots", "02_quality_control")
processed_dir = here("processed-data", "02_quality_control")

location <- rowRanges(sce)
is.mito <- any(rownames(sce) == "MT")

library(scuttle)
df <- perCellQCMetrics(sce, subsets=list(Mt=grep("^MT-", rownames(sce))))
summary(df$subsets_Mito_sum)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 503    3124    5039    5776    7696   38904

sce <- addPerCellQCMetrics(sce, subsets=list(Mito=is.mito))
colnames(colData(sce))
# [1] "sum"                   "detected"              "subsets_Mito_sum"      "subsets_Mito_detected"
# [5] "subsets_Mito_percent"  "total"

sum(colData(sce)$subsets_Mito_sum)
# [1] 0 

# reasons <- perCellQCFilters(df,
#                             sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))
# colSums(as.matrix(reasons))

# ===== Problem: I'm getting zero mitochondrial reads =====
# let's test this by searching for mito genes without the 'MT-'

mitochondrial_genes_list <- c("ND1", "ND2", "COX1", "COX2-201", "ATP6", "CYTB", "ATP6")

existing_mito_genes <- mitochondrial_genes_list[mitochondrial_genes_list %in% rownames(sce)]
print(existing_mito_genes)
# character(0)

pdf(width=10, height=5, here(plot_dir, "UMI_violin_macaque.pdf"))
plotColData(sce, y="detected", x="Sample", color_by='Sample') +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

mean(colSums(counts(sce)))
# [1] 7983.276

median(colSums(counts(sce)))
# [1] 3819


# plot he number of cells per sample
sample_metadata <- colData(sce)$Sample
cells_per_sample <- table(sample_metadata)
plot_data <- as.data.frame(cells_per_sample)

ggplot(plot_data, aes(x=sample_metadata, y=Freq)) +  # Var1 contains sample names, Freq contains cell counts
  geom_bar(stat="identity", fill="blue", alpha=0.7) +
  labs(title="Number of Cells per Sample",
       x="Sample",
       y="Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Optional: Rotate x-axis labels for better readability

ggsave(here(plot_dir, "Nuclei_per_sample_macaque.pdf"))


# ====== Normalization ======

clust.sce <- quickCluster(sce) 
sce <- computeSumFactors(sce, cluster=clust.sce, min.mean=0.1)
sce <- logNormCounts(sce)



# =========== Dim Reduction ============

# ===== PCA =====
library(scran)
top.genes <- getTopHVGs(sce, n=2000)

set.seed(100) # See below.
sce <- fixedPCA(sce, subset.row=top.genes) 
reducedDimNames(sce)


# plot by sample info
features = c("Region", "Batch", "Anatomoical")

pdf(here(plot_dir,"PCA_macaque_by_sampleinfo.pdf"))
for (feature in features) {
p <- plotPCA(sce, ncomponents=5, colour_by=feature)
print(p)
}
dev.off()


# plot by genes
features = c("SNAP25", "MBP", "GAD1", "SLC17A7","SLC1A2", "GULP1", 'COL25A1','CRH')

pdf(here(plot_dir,"PCA_macaque_by_genes.pdf"))
for (feature in features) {
  p <- plotPCA(sce, colour_by=feature)
  print(p)
}
dev.off()


# ===== TSNE =====
sce <- runTSNE(sce, dimred="PCA")

# plot by sample info 
features = c("Sample", "Region", "Batch", "Anatomoical")

pdf(here(plot_dir,"TSNE_macaque_by_sampleinfo.pdf"))
for (feature in features) {
  p <- plotReducedDim(sce, dimred="TSNE", color_by=feature)
  print(p)
}
dev.off()

# plot by genes
features = c("SNAP25", "MBP", "GAD1", "SLC17A7","SLC1A2", "GULP1", 'COL25A1','CRH')

pdf(here(plot_dir,"TSNE_macaque_by_genes.pdf"))
for (feature in features) {
  p <- plotReducedDim(sce, dimred="TSNE", color_by=feature)
  print(p)
}
dev.off()

# ===== UMAP =====
sce <- runUMAP(sce, dimred="PCA")

# plot by sample info 
features = c("Sample", "Region", "Batch", "Anatomoical")

pdf(here(plot_dir,"UMAP_macaque_by_sampleinfo.pdf"))
for (feature in features) {
  p <- plotReducedDim(sce, dimred="UMAP", color_by=feature)
  print(p)
}
dev.off()

# plot by genes
features = c("SNAP25", "MBP", "GAD1", "SLC17A7","SLC1A2", "GULP1", 'COL25A1','CRH')

pdf(here(plot_dir,"UMAP_macaque_by_genes.pdf"))
for (feature in features) {
  p <- plotReducedDim(sce, dimred="UMAP", color_by=feature)
  print(p)
}
dev.off()


# ========== Harmony =========

# === Sample correction ===
sce.harmony.sample <- RunHarmony(sce, group.by.vars = "Sample", verbose = TRUE)

# Run TSNE on Harmony.Sample
sce.harmony.sample <- runTSNE(sce.harmony.sample, dimred = "HARMONY")

# Plot TSNE x Sample info with harmony corrected for sample
pdf(here(plot_dir,"TSNE_macaque_by_sampleinfo_harmony_bysample.pdf"))
for (feature in features) {
  p <- plotReducedDim(sce.harmony.sample, dimred="TSNE", color_by=feature)
  print(p)
}
dev.off()

# === Batch correction ===
sce.harmony.batch <- RunHarmony(sce, group.by.vars = "Batch", verbose = TRUE)

# Run TSNE on Harmony.Sample
sce.harmony.batch <- runTSNE(sce.harmony.batch, dimred = "HARMONY")

# Plot TSNE x Sample info with harmony corrected for sample
pdf(here(plot_dir,"TSNE_macaque_by_sampleinfo_harmony_bybatch.pdf"))
for (feature in features) {
  p <- plotReducedDim(sce.harmony.batch, dimred="TSNE", color_by=feature)
  print(p)
}
dev.off()




# ====== clustering =======

message("running buildSNNGraph - ", Sys.time())
snn.gr.25 <- buildSNNGraph(sce, k = 25, use.dimred = "PCA")
snn.gr.50 <- buildSNNGraph(sce, k = 50, use.dimred = "PCA")

message("running walktrap - ", Sys.time())
clust20 <- igraph::cluster_walktrap(snn.gr.25)$membership
clust50 <- igraph::cluster_walktrap(snn.gr.50)$membership
