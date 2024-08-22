# originally copied from https://github.com/LieberInstitute/spatialdACC/blob/main/code/17_LDSC/spatial/01_aggregate.R

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)

sce <- readRDS(here("processed-data", "09_psychiatric_GSEA", "sce.human_all_genes.rds"))

LA_clusters <- c("GULP1_TRHDE", "ZBTB20_SLC4A4")
BA_clusters <- c("PEX5L_MYRIP", "MEIS2_COL25A1")
aBA_clusters <- c("ESR1_ADRA1A", "GRIK3_TNS3")

sce$subregion_celltype <- "Other"
sce$subregion_celltype[sce$fine_celltype %in% LA_clusters] <- "LA"
sce$subregion_celltype[sce$fine_celltype %in% BA_clusters] <- "BA"
sce$subregion_celltype[sce$fine_celltype %in% aBA_clusters] <- "aBA"

sce_pseudo <- aggregateAcrossCells(sce, id=colData(sce)[,c("subregion_celltype")])
sce_pseudo <- sce_pseudo[,sce_pseudo$subregion_celltype != "Other"]
dim(sce_pseudo)


#find a good expression cutoff using edgeR::filterByExpr
rowData(sce_pseudo)$high_expr_group_cluster <- filterByExpr(sce_pseudo, group = sce_pseudo$subregion_celltype)

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
colnames(normd)<-sce_pseudo$subregion_celltype

write.table(normd, here::here("processed-data", "09_psychiatric_GSEA","LDSC", "specificity_score_putative", "human_sce_aggregated.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
