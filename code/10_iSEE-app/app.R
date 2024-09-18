library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")
library("scuttle")
library("SpatialExperiment")

#Load the object
load("sce_clean.rda", verbose = TRUE)

#Change the rownames frome ensembl id to gene_name
rownames(sce) <-  uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)

#Source
source("initial.R", print.eval = TRUE)

#Deploy app
iSEE(
  sce,
  appTitle = "LS snRNAseq data",
  initial = initial
)
