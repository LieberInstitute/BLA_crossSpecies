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
processed_dir = here("processed-data", "05_species_comparisons")
plot_dir = here("plots", "05_species_comparisons")

Amy.700 <- readRDS(here("processed-data", "00_costa_rds", "amygdalaCNScaled700.rds"))
Amy.700
sce.mac <- as.SingleCellExperiment(Amy.700, assay="RNA")
sce.mac
# class: SingleCellExperiment 
# dim: 19353 120794 
# metadata(0):
#   assays(2): counts logcounts
# rownames(19353): ZNF692 ZNF672 ... LILRA6 USP29
# rowData names(0):
#   colnames(120794): BasalVentralAP2_AAACCCAGTAGGTCAG-1 BasalVentralAP2_AAACCCAGTCAAATCC-1 ...
# Central_Nucleus_AP3_TTTGGTTGTGGCTAGA-1 Central_Nucleus_AP3_TTTGTTGCACATGGTT-1
# colData names(11): orig.ident nCount_RNA ... seurat_clusters ident
# reducedDimNames(2): PCA UMAP
# mainExpName: RNA
# altExpNames(0):


Bab.200 <-readRDS(here("processed-data", "00_costa_rds", "AllBaboon_Amygdala200_04.rds"))
Bab.200
sce.bab <- as.SingleCellExperiment(Bab.200, assay="RNA")
sce.bab
# class: SingleCellExperiment 
# dim: 29544 71903 
# metadata(0):
#   assays(2): counts logcounts
# rownames(29544): TMEM88B ANKRD65 ... LOC101007301 LOC103884182
# rowData names(0):
#   colnames(71903): Baboon2_Basal1_AAACCCAAGACTGTTC-1 Baboon2_Basal1_AAACCCAAGATCCAAA-1 ... An3_LA_L2_TTTGTTGAGTGCTAGG-1
# An3_LA_L2_TTTGTTGTCTCCGATC-1
# colData names(8): orig.ident nCount_RNA ... seurat_clusters ident
# reducedDimNames(2): PCA UMAP
# mainExpName: RNA
# altExpNames(0):

load(here("processed-data", "00_costa_rds", "sce_annotated.rda"))
sce.totty <- sce
# class: SingleCellExperiment 
# dim: 36601 19682 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(36601): MALAT1 ERBB4 ... AC136616.2 AC141272.1
# rowData names(8): source type ... Symbol.uniq binomial_deviance
# colnames(19682): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ... 5_TTTGTTGTCGGACTGC-1 5_TTTGTTGTCGTTGTTT-1
# colData names(24): Sample Barcode ... celltype annotation
# reducedDimNames(4): PCA TSNE UMAP HARMONY
# mainExpName: NULL
# altExpNames(0):


# ======== Normalize sce.totty ==========
#apparently sce.totty might not be normalized. Let's try to convert to seuratr and normalize the same way Vinny did
any_negatives_human <- any(logcounts(sce.totty) < 0)

seurat.totty <- as.Seurat(sce.totty, counts = "counts", data = NULL)
seurat.totty <- NormalizeData(object = seurat.totty)
sce.totty <- as.SingleCellExperiment(seurat.totty)


logcounts(sce.totty)
logcounts(sce.mac)


# ====== Annotate NHP ======

# baboon
sce.bab$broad_celltype[sce.bab$ident %in% c(14,20,25,18,13,10,15)]<-'Inhib'
sce.bab$broad_celltype[sce.bab$ident %in% c(2,1,0,22,12,4,8,17,19,9,26,11,16)]<-'Excit'
sce.bab$broad_celltype[sce.bab$ident %in% c(5)]<-'Micro'
sce.bab$broad_celltype[sce.bab$ident %in% c(27)]<-'Endo'
sce.bab$broad_celltype[sce.bab$ident %in% c(23,3)]<-'Oligo'
sce.bab$broad_celltype[sce.bab$ident %in% c(6,24)]<-'Astro'
sce.bab$broad_celltype[sce.bab$ident %in% c(21,7)]<-'OPC'

# macaque
sce.mac$broad_celltype[sce.mac$ident %in% c(12,10,22,31,30,15,24,19,27)]<-'Inhib'
sce.mac$broad_celltype[sce.mac$ident %in% c(23,7,13,9,5,8,4,16,14,20,28,11)]<-'Excit'
sce.mac$broad_celltype[sce.mac$ident %in% c(6)]<-'Micro'
sce.mac$broad_celltype[sce.mac$ident %in% c(21,17)]<-'Endo'
sce.mac$broad_celltype[sce.mac$ident %in% c(18,0,2)]<-'Oligo'
sce.mac$broad_celltype[sce.mac$ident %in% c(3,26)]<-'Astro'
sce.mac$broad_celltype[sce.mac$ident %in% c(1,25,29)]<-'OPC'

# ====== Make broader celltypes for Totty data ======
celltypes <- colData(sce.totty)$annotation

# Update based on the specific keywords
broad_celltypes <- ifelse(grepl("Exc", celltypes), "Excit", 
                          ifelse(grepl("Inh", celltypes), "Inhib", 
                                 ifelse(grepl("Astro", celltypes), "Astro", celltypes)))

# Update the broad_celltype column in the SCE object
colData(sce.totty)$broad_celltype <- broad_celltypes



# ================ Combine datasets ================

# Add species information to each SCE
colData(sce.totty)$species <- "human"
colData(sce.bab)$species <- "baboon"
colData(sce.mac)$species <- "macaque"


# drop uneven assays
#assays(sce.totty)$counts <- NULL
assays(sce.totty)$binomial_deviance_residuals <- NULL
rowData(sce.totty) <- NULL
metadata(sce.totty) <- list()



# drop all reduced dims
reducedDims(sce.totty) <- NULL
reducedDims(sce.mac) <- NULL
reducedDims(sce.bab) <- NULL

# drop alternate experiment names
altExp(sce.bab) <- NULL
altExp(sce.mac) <- NULL

# Get the shared genes 
shared_genes <- Reduce(intersect, list(rownames(sce.totty), rownames(sce.bab), rownames(sce.mac)))

sce.totty <- sce.totty[shared_genes, ]
sce.bab <- sce.bab[shared_genes, ]
sce.mac <- sce.mac[shared_genes, ]


# get the shared columns
shared_colnames <- Reduce(intersect, list(colnames(colData(sce.totty)), 
                                          colnames(colData(sce.bab)), 
                                          colnames(colData(sce.mac))))

colData(sce.totty) <- colData(sce.totty)[, shared_colnames]
colData(sce.bab) <- colData(sce.bab)[, shared_colnames]
colData(sce.mac) <- colData(sce.mac)[, shared_colnames]


# combine
combined_counts <- cbind(assay(sce.mac, "counts"), 
                            assay(sce.bab, "counts"), 
                            assay(sce.totty, "counts"))

combined_colData <- rbind(colData(sce.mac), 
                          colData(sce.bab), 
                          colData(sce.totty))


combined_sce <- SingleCellExperiment(assays=list(counts=combined_counts), colData=combined_colData)
saveRDS(combined_sce, file = here(processed_dir, "combined_sce.rds"))


# =========== Normalization
combined_sce <- readRDS(here('processed-data', "05_species_comparisons", "combined_sce.rds"))



combined_sce <- logNormCounts(combined_sce)


# ========= Feature (highly variable gene) selection ==========

dec.pbmc <- modelGeneVar(combined_sce)

# Visualizing the fit:
fit.pbmc <- metadata(dec.pbmc)
plot(fit.pbmc$mean, fit.pbmc$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.pbmc$trend(x), col="dodgerblue", add=TRUE, lwd=2)

# Ordering by most interesting genes for inspection.
dec.pbmc[order(dec.pbmc$bio, decreasing=TRUE),] 

# get top 2000 HVGs
chosen <- getTopHVGs(dec.pbmc, n=5000)

combined_sce.hvg <- combined_sce[chosen,]
dim(combined_sce.hvg)


# ====== Pseudobulking =======
# aggregate both the raw counts and logcounts

combined_sce.psb <- aggregateAcrossCells(combined_sce.hvg, id=colData(combined_sce.hvg)[,c("broad_celltype", "species")], use.assay.type=c("counts","logcounts"))
combined_sce.psb



#======== Before running PCA ...... ============
#
# Thus far, the PCA1-2 was dominated by variance in NHP excitatory celltypes,
# and macaque Oligos.
#
# Stephanie pointed out that this is likely due to library size. We need to normalize this. 
#
#
# To first get a sense for this, lets plot MBP and SNAP25 to get a sense of Oligo and Excit library sizes
#
#

pdf(here(plot_dir,"Expression_pseudobulk_MBP.pdf"))
plotExpression(combined_sce.psb, features="MBP", x='species', colour_by='broad_celltype')
dev.off()

pdf(here(plot_dir,"Expression_pseudobulk_SNAP25.pdf"))
plotExpression(combined_sce.psb, features="SNAP25", x='species', colour_by='broad_celltype')
dev.off()

# Results:
# Yep, these exact same cell types show a hugely inflated number of counts. 
# we need to normalize this. somehow.
#

# ========= Normalizing with size factors using DESeq2 ==========
#
# We will first need to convert our SCE to a DEseq object.
#


library(DESeq2)
counts <- counts(combined_sce.psb)
sampleInfo <- colData(combined_sce.psb)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleInfo, design = ~ c(species, broad_celltype))

dds <- estimateSizeFactors(dds)

# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))
dds_norm <- DESeqTransform(se)

colnames(dds_norm) <- c(1:21)

pdf(here(plot_dir,"DESeq2_PCA.pdf"))
pcaData <- DESeq2::plotPCA(dds_norm, ntop = 500, intgroup = c("species", "broad_celltype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=broad_celltype, shape=species)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

# Plot PCA
pdf(here(plot_dir,"DESeq2_PCA_rlog.pdf"))
pcaData <- DESeq2::plotPCA(rld, ntop = 500, intgroup = c("species", "broad_celltype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=broad_celltype, shape=species)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()



# ===== Heatmaps for hierarchical clustering ======

# Extract sample-level variables
metadata <- colData(combined_sce.psb) %>% 
  as.data.frame() %>% 
  dplyr::select(species, broad_celltype)
metadata
# species broad_celltype
# 1   baboon          Astro
# 2    human          Astro
# 3  macaque          Astro
# 4   baboon           Endo
# 5    human           Endo
# 6  macaque           Endo
# 7   baboon          Excit
# 8    human          Excit
# 9  macaque          Excit
# 10  baboon          Inhib
# 11   human          Inhib
# 12 macaque          Inhib
# 13  baboon          Micro
# 14   human          Micro
# 15 macaque          Micro
# 16  baboon          Oligo
# 17   human          Oligo
# 18 macaque          Oligo
# 19  baboon            OPC
# 20   human            OPC
# 21 macaque            OPC

# Define custom colors for broad_celltype
broad_celltype_colors <- c(
  "Oligo" = "#E69F00",
  "Inhib" = "#56B4E9",
  "Excit" = "#009E73",
  "Astro" = "#F0E442",
  "Micro" = "#0072B2",
  "OPC" = "#D55E00",
  "Endo" = "#CC79A7"
)

# Define custom colors for species
species_colors <- c(
  "macaque" = "#D7263D",
  "human" = "#26D5B8",
  "baboon" = "#F46036"
)

annotation_colors <- list(
  broad_celltype = broad_celltype_colors,
  species = species_colors
)

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap
pdf(here(plot_dir,"DESeq2_heatmap_rlog.pdf"))
pheatmap(rld_cor, annotation=metadata, annotation_colors=annotation_colors)
dev.off()







