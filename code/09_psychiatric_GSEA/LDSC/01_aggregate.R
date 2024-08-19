# originally copied from https://github.com/LieberInstitute/spatialdACC/blob/main/code/17_LDSC/spatial/01_aggregate.R

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)

sce <- readRDS(here("processed-data", "09_psychiatric_GSEA", "sce.human_all_genes.rds"))

# filter cell types with < 70 cells
celltype_counts <- table(colData(sce)$fine_celltype)
celltypes_to_drop <- names(celltype_counts[celltype_counts < 70])
sce <- sce[, !(colData(sce)$fine_celltype %in% celltypes_to_drop)]
sce$fine_celltype <- factor(sce$fine_celltype)

sce_pseudo <- aggregateAcrossCells(sce, id=colData(sce)[,c("fine_celltype")])
dim(sce_pseudo)

#find a good expression cutoff using edgeR::filterByExpr
rowData(sce_pseudo)$high_expr_group_cluster <- filterByExpr(sce_pseudo, group = sce_pseudo$fine_celltype)

summary(rowData(sce_pseudo)$high_expr_group_cluster)
#   Mode   FALSE    TRUE 
# logical     600   21244 

## Now filter
sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$high_expr_group_cluster, ]
sce_pseudo <- sce_pseudo[!duplicated(rowData(sce_pseudo)$gene_name),]
sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$gene_type=='protein_coding',]
dim(sce_pseudo)
# [1] 15128    30

normd <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo))
rownames(normd)<-rowData(sce_pseudo)$gene_name
colnames(normd)<-sce_pseudo$fine_celltype

write.table(normd, here::here("processed-data", "09_psychiatric_GSEA","LDSC", "human_sce_aggregated.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
