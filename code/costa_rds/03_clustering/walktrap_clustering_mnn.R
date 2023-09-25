# qrsh -l mem_free=120G,h_vmem=120G -now n

library("SingleCellExperiment")
#library("jaffelab")
library("scater")
library("scran")
library("here")
library("sessioninfo")

# Save directories
plot_dir = here("plots", "costa_rds", "03_clustering")
processed_dir = here("processed-data","costa_rds", "03_clustering")

load(here("processed-data", "costa_rds",  "03_clustering", "sce_combined_harmony.rda"), verbose = TRUE)
load(here("processed-data",  "costa_rds", "03_clustering", "sce_combined_mnn.rda"), verbose = TRUE)
sce <- mnn.out

# Extract the "corrected" reduced dimension from mnn.out
corrected_dim <- reducedDim(mnn.out, "corrected")

# Add the extracted dimension to the sce object
reducedDim(sce, "MNN") <- corrected_dim
rm(mnn.out)

# ============== Buld SNN Graph ============
#
# READ ME: I added the MNN reduced dim to sce from the Harmony output
# On second thought, this might not make an sense at all. Does Harmony/MNN
# change the underlying cout/logcount values?
#
message("running buildSNNGraph - ", Sys.time())
snn.gr.50 <- buildSNNGraph(sce, k = 50, use.dimred = "MNN")

message("running walktrap - ", Sys.time())
clust50 <- igraph::cluster_walktrap(snn.gr.50)$membership

table(clust50)


##add to sce
sce$k_50_label<-clust50

message("saving data - ", Sys.time())
save(sce, file=here(processed_dir, "sce_clustered_k50.rda"))

##make some prelim plots
pdf(here(plot_dir,"UMAP_k50_mnn.pdf"))
plotUMAP(sce,colour_by='k_50_label',text_by='k_50_label')
dev.off()

##make some prelim plots
pdf(here(plot_dir,"TSNE_k50_mnn.pdf"))
plotTSNE(sce,colour_by='k_50_label',text_by='k_50_label')
dev.off()



##add logcounts for viz
message("normalizing counts - ", Sys.time())
set.seed(1000)
sce <- computeSumFactors(sce, cluster=sce$k_50_label)
sce <- logNormCounts(sce)
sce

message("saving data - ", Sys.time())
save(sce, file=here(processed_dir, "sce_clustered_k50_normalized.rda"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


