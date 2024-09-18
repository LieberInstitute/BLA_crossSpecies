library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")
library("scuttle")
library("SpatialExperiment")

#Load the object
load("sce_FINAL_all_celltypes.rda", verbose = TRUE)
#sce <- rda

#Change the rownames frome ensembl id to gene_name
#rownames(sce) <-  uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)

#Source
source("initial.R", print.eval = TRUE)

#Deploy app
iSEE(
  rda,
  appTitle = "BLA CrossSpecies",
  initial = initial
)


