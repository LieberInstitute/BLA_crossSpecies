# qrsh -l mem_free=80G,h_vmem=80G -now n
library("SingleCellExperiment")
#library("jaffelab")
library("scater")
library("scran")
library("here")
library("sessioninfo")

# Save directories
plot_dir = here("plots", "03_clustering")
processed_dir = here("processed-data","03_clustering")

load(here("processed-data", "03_clustering", "sce_combined_harmony.rda"), verbose = TRUE)
load(here("processed-data", "03_clustering", "sce_combined_mnn.rda"), verbose = TRUE)
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

##make some prelim plots
pdf(here(plot_dir,"UMAP_k50_mnn.pdf"))
plotUMAP(sce,colour_by='k_10_label',text_by='k_10_label')
dev.off()

##make some prelim plots
pdf(here(plot_dir,"TSNE_k50_mnn.pdf"))
plotTSNE(sce,colour_by='k_10_label',text_by='k_10_label')
dev.off()



##add logcounts for viz
message("normalizing counts - ", Sys.time())
set.seed(1000)
sce <- computeSumFactors(sce, cluster=sce$k_50_label)
sce <- logNormCounts(sce)
sce

message("saving data - ", Sys.time())
save(sce, file=here(processed_dir, "sce_clustered_Ks_10_60.rda"))

pdf(here(plot_dir,"markers_k50.pdf"))
plotExpression(sce,features=c('SYT1','SLC17A7','SLC17A6',
                              'GAD1','GAD2','MBP',
                              'GFAP','TTR','CSF1R',
                              'PDGFRA','FLT1','TNNT2'),
               x="k_50_label", colour_by="k_50_label", point_alpha=0.5, point_size=.7,add_legend=F)+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.3)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

