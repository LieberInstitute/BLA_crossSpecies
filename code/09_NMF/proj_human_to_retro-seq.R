library(SpatialExperiment)
library(scater)
library(RcppML)
library(ggspavis)
library(here)
library(scRNAseq)
library(Matrix)
library(scran)
library(scuttle)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(bluster)
library(patchwork)
library(cowplot)
library(projectR)
library(spatialLIBD)

plot_dir <- here("plots","09_NMF")
processed_dir <- here("processed-data")

# load NMF results
load(here(processed_dir,"09_NMF", "RcppML_NMF_LIBD.rda"))
x

# load retro-seq object
sce.retro <- readRDS(here("processed-data","retro-seq-data", "retro-seq_sce_final.rds"))
sce.retro
# class: SingleCellExperiment 
# dim: 19507 1466 
# metadata(1): log
# assays(1): logcounts
# rownames(19507): Xkr4 Rp1 ... mt-Nd6 mt-Cytb
# rowData names(3): chrom gene_name human_ortholog
# colnames(1466): Pool154_Plate1-1-A1-A13 Pool154_Plate1-1-A1-A14 ...
# Pool180_Plate13-6-I24-P23 Pool180_Plate13-6-I24-P24
# colData names(21): cell.name Cell.Name ... Subclass graph_clust
# reducedDimNames(3): PCA UMAP HARMONY
# mainExpName: NULL
# altExpNames(0):

# load LIBD data
load(here("processed-data","retro-seq-data", "BLA_sce_annotated.rda"))
sce.libd <- sce
sce.libd
# class: SingleCellExperiment 
# dim: 36601 21268 
# metadata(1): Samples
# assays(3): counts binomial_deviance_residuals logcounts
# rownames(36601): MALAT1 ERBB4 ... AC136616.2 AC141272.1
# rowData names(8): source type ... Symbol.uniq binomial_deviance
# colnames(21268): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ... 5_TTTGTTGTCGGACTGC-1
# 5_TTTGTTGTCGTTGTTT-1
# colData names(23): Sample Barcode ... k_50_label celltype
# reducedDimNames(4): PCA TSNE UMAP HARMONY
# mainExpName: NULL
# altExpNames(0):




# ====== Set up patterns and loadings ========

# extract patterns
patterns <- t(x$h)
colnames(patterns) <- paste("NMF", 1:100, sep = "_")

loadings <- x$w
rownames(loadings) <- rownames(sce.libd)


# ====== Subset to only 1:1 orthologs ======
library(orthogene)
method <- "gprofiler"

mouse_orthos <- orthogene::convert_orthologs(gene_df = logcounts(sce.retro),
                                        gene_input = "rownames",
                                        gene_output = "rownames",
                                        input_species = "mouse",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species",
                                        method = method)
# =========== REPORT SUMMARY ===========
#     
# Total genes dropped after convert_orthologs :
#     3,242 / 19,493 (17%)
# Total genes remaining after convert_orthologs :
#     16,251 / 19,493 (83%)

# subset to remaining genes and swap to human names
sce.retro <- sce.retro[rownames(sce.retro) %in% names(rownames(mouse_orthos))]
rowData(sce.retro)$gene_name <- rownames(sce.retro)
rowData(sce.retro)$human_ortholgs <- rownames(mouse_orthos)
rownames(sce.retro) <- rowData(sce.retro)$human_ortholgs



# get only orthologs that are in the human dataset
valid.human.orthologs <- rownames(mouse_orthos)[rownames(mouse_orthos) %in% rownames(sce.libd)]
length(valid.human.orthologs)
# [1] 15938

# subset to only 1:1 orthologs
sce.libd <- sce.libd[valid.human.orthologs ,]
sce.libd
# class: SingleCellExperiment
# dim: 13874 21268
# metadata(1): Samples
# assays(1): counts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(7): source type ... gene_type Symbol.uniq
# colnames(21268): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 5_TTTGTTGTCGGACTGC-1 5_TTTGTTGTCGTTGTTT-1
# colData names(27): Sample Barcode ... discard_auto doubletScore
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):


# ====== project loadings to retro-seq data =======
# drop any rownames in SPE not in scee
sce.retro <- sce.retro[rownames(sce.retro) %in% rownames(sce.libd),]
loadings <- loadings[rownames(loadings) %in% rownames(sce.libd),]

# drop any rownames in loadings not in spe
loadings <- loadings[rownames(loadings) %in% rownames(sce.retro),]
sce.retro <- sce.retro[rownames(sce.retro) %in% rownames(loadings),]
sce.retro <- sce.retro[!duplicated(rownames(sce.retro)),]

logcounts <- logcounts(sce.retro)
logcounts <- Matrix(logcounts, sparse = TRUE)

proj <- project(logcounts, loadings)
proj <- t(proj)
colnames(proj) <- paste("NMF", 1:100, sep = "_")



# add to reducedDims
reducedDim(sce.retro, "NMF_proj") <- proj

spe.temp <- sce.retro

# add each proj column to colData(spe)
for (i in 1:100){
    colData(spe.temp)[[paste0("NMF_",i)]] <- reducedDims(spe.temp)$NMF_proj[,i]
}

# ======== Heatmaps ========

# ==== Spatial domains ====
# drop NA subclass

# create dataframe 
data <- data.frame(colData(sce.retro), reducedDims(sce.retro)$NMF_proj)

# aggregate NMF patterns across clusters. # grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("NMF", colnames(data))],
                      by=list(data$Cocluster), 
                      FUN=mean)

# move Group.1 to row names, then drop
rownames(agg_data) <- agg_data$Group.1
agg_data <- agg_data[,-1]

p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("grey","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               #scale="column"
)

pdf(here(plot_dir, "Heatmap_NMF_retro-seq.pdf"), width=12.5, height=5)
p1
dev.off()



# ======= Exporting top gene per factor ========

# function for getting top n genes for each pattern
top_genes <- function(W, n=10){
    top_genes <- apply(W, 2, function(x) names(sort(x, decreasing=TRUE)[1:n]))
    return(top_genes)
}

#get top genes
top15 <- top_genes(loadings, 15)
head(top15)

# add colnames
colnames(top15) <- colnames(patterns)

# export to csv
write.csv(top15, here(processed_dir, "NMF_Yu", "top15_genes_per_factor.csv"))
