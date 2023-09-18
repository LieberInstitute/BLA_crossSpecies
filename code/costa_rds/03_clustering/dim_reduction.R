library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("here")
library("sessioninfo")
library("harmony")
library("BiocSingular")

# Save directories
plot_dir = here("plots", "03_clustering")
processed_dir = here("processed-data", "03_clustering")

## Load sce
sce <- readRDS(here('processed-data', "05_species_comparisons", "combined_sce.rds"))
sce <- logNormCounts(sce)

# ========= Feature (highly variable gene) selection ==========

dec.pbmc <- modelGeneVar(sce)

# Visualizing the fit:
fit.pbmc <- metadata(dec.pbmc)
plot(fit.pbmc$mean, fit.pbmc$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.pbmc$trend(x), col="dodgerblue", add=TRUE, lwd=2)

# Ordering by most interesting genes for inspection.
dec.pbmc[order(dec.pbmc$bio, decreasing=TRUE),] 

# get top 2000 HVGs
chosen <- getTopHVGs(dec.pbmc, n=2000)

# =========== PCA ===========

set.seed(915)
message("running PCA - ", Sys.time())
sce <- fixedPCA(sce, subset.row=chosen)


# =========== TASNE ===========

message("running TSNE - ", Sys.time())
sce <- runTSNE(sce, dimred = "PCA")

pdf(here(plot_dir,"TSNE_combined_uncorrected.pdf"))
plotTSNE(sce, colour_by='species')
dev.off()

set.seed(916)
message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "PCA")

pdf(here(plot_dir,"UMAP_combined_uncorrected.pdf"))
plotUMAP(sce,colour_by='Sample')
dev.off()

message("Saving Data - ", Sys.time())
save(sce, file = here(processed_dir, "sce.rda"))
# sgejobs::job_single('normalize_step1_glm', create_shell = TRUE, queue= 'bluejay', memory = '75G', command = "Rscript normalize_step1_glm.R")

# Run harmony
message("running Harmony - ", Sys.time())
sce <- RunHarmony(sce, group.by.vars = "species", verbose = TRUE)


#### TSNE & UMAP ####
set.seed(602)
message("running TSNE - ", Sys.time())
sce <- runTSNE(sce, dimred = "HARMONY")

pdf(here(plot_dir,"TSNE_HarmonyCorrected_by_species.pdf"))
plotTSNE(sce,colour_by='species')
dev.off()

message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "HARMONY")

pdf(here(plot_dir,"UMAP_HarmonyCorrected_by_species.pdf"))
plotUMAP(sce,colour_by='species')
dev.off()

message("Done UMAP - Saving data...", Sys.time())
save(sce, file = here(processed_dir, "sce_combined_harmony.rda"))



# ============== Batch correction with MNN ================
library(batchelor)

set.seed(1000101001)
mnn.out <- fastMNN(sce, batch=sce$species, d=50, k=20, subset.row=chosen,
                   BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
mnn.out

#### TSNE & UMAP ####
set.seed(602)
message("running TSNE - ", Sys.time())
mnn.out <- runTSNE(mnn.out, dimred = "corrected")

pdf(here(plot_dir,"TSNE_corrected_by_species.pdf"))
plotTSNE(mnn.out,colour_by='batch')
dev.off()

message("running UMAP - ", Sys.time())
mnn.out <- runUMAP(mnn.out, dimred = "corrected")

pdf(here(plot_dir,"UMAP_corrected_by_species.pdf"))
plotUMAP(mnn.out,colour_by='batch')
dev.off()

message("Done UMAP - Saving data...", Sys.time())
save(mnn.out, file = here(processed_dir, "sce_combined_mnn.rda"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

