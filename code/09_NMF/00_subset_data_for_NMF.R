library(SingleCellExperiment)  
library(here)                

# ====== Load Retro seq and Excitatory data ======
sce.retro <- readRDS(here("processed-data","retro-seq-data", "retro-seq_sce_final.rds"))
sce.retro
# class: SingleCellExperiment 
# dim: 19507 1466 
# metadata(1): log
# assays(1): logcounts
# rownames(19507): Xkr4 Rp1 ... mt-Nd6 mt-Cytb
# rowData names(3): chrom gene_name human_ortholog
# colnames(1466): Pool154_Plate1-1-A1-A13 Pool154_Plate1-1-A1-A14 ... Pool180_Plate13-6-I24-P23
# Pool180_Plate13-6-I24-P24
# colData names(21): cell.name Cell.Name ... Subclass graph_clust
# reducedDimNames(3): PCA UMAP HARMONY
# mainExpName: NULL
# altExpNames(0):


sce.excit <- readRDS(here("processed-data", "07_annotation", "sce.excit.integrated.annotated.rds"))
sce.excit
# class: SingleCellExperiment 
# dim: 13874 70227 
# metadata(0):
#     assays(2): counts logcounts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(0):
#     colnames(70227): AAACCCACAGGTCCCA-1 AAACCCATCCATTGTT-1 ... TGGGAAGGTTAGCGGA-1 TTCACGCGTAGTCTTG-1
# colData names(61): orig.ident nCount_originalexp ... k.60_cluster.fun.louvain fine_type
# reducedDimNames(2): PCA UMAP
# mainExpName: originalexp
# altExpNames(0):

# ====== Subset data for NMF ======
# We need to get only excitatory neuron types from mouse retro data.
# Also, we need to subset to just human samples from sce.excit

# == Retro seq data ==
sce.retro$Subclass <-sce.retro$Subclass$Subclass # was a Data Frame with two columns for some reason

# subset to only Subclass that contains "Glut"
sce.retro.excit <- sce.retro[, grepl("Glut", sce.retro$Subclass)]
sce.retro.excit
# class: SingleCellExperiment 
# dim: 19507 1137 
# metadata(1): log
# assays(1): logcounts
# rownames(19507): Xkr4 Rp1 ... mt-Nd6 mt-Cytb
# rowData names(3): chrom gene_name human_ortholog
# colnames(1137): Pool154_Plate1-1-A1-A13 Pool154_Plate1-1-A1-A14 ... Pool180_Plate13-6-I24-P23
# Pool180_Plate13-6-I24-P24
# colData names(21): cell.name Cell.Name ... Subclass graph_clust
# reducedDimNames(3): PCA UMAP HARMONY
# mainExpName: NULL
# altExpNames(0):

# save
saveRDS(sce.retro.excit, file = here("processed-data", "09_NMF", "sce.retro.excit.rds"))

# == Human data ==
sce.excit.human <- sce.excit[, sce.excit$species == "human"]
sce.excit.human
# class: SingleCellExperiment 
# dim: 13874 12160 
# metadata(0):
#     assays(2): counts logcounts
# rownames(13874): ANKRD65 AURKAIP1 ... R3HDM4 KISS1R
# rowData names(0):
#     colnames(12160): AAACCCACAGGTCCCA-1 AAACCCATCCATTGTT-1 ... TTTGGTTAGCTCGACC-1 TTTGTTGTCGTTGTTT-1
# colData names(61): orig.ident nCount_originalexp ... k.60_cluster.fun.louvain fine_type
# reducedDimNames(2): PCA UMAP
# mainExpName: originalexp
# altExpNames(0):

# save
saveRDS(sce.excit.human, file = here("processed-data", "09_NMF", "sce.excit.human.rds"))
