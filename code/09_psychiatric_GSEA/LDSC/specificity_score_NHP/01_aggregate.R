# originally copied from https://github.com/LieberInstitute/spatialdACC/blob/main/code/17_LDSC/spatial/01_aggregate.R

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)


sce.excit <- readRDS(here("processed-data", "07_annotation_and_characterization", "sce_excit_final_subclusters_annotated.rds"))

# ======== Pseudobulk analysis ========
# for DE analysis, it makes the most sense to subset to anatomical samples
# and perform pseudobulk and DE analysis. This avoids the statistical "double dipping"
# of using markers genes based on clusters... which are based on DEGs.

# subset SCE to just Later and Basal
sce.subset <- sce.excit[, sce.excit$Subregion %in% c("Lateral", "Basal", "Accessory Basal")]
dim(sce.subset)
# [1] 13874 53306

# aggregate across cells
sce_pseudo <- aggregateAcrossCells(sce.subset, ids = colData(sce.subset)[,c("Subregion")])
dim(sce_pseudo)
#[1] 13874    3

unique(sce_pseudo$species)
# [1] "baboon"  "macaque"

#find a good expression cutoff using edgeR::filterByExpr
rowData(sce_pseudo)$high_expr_group_cluster <- filterByExpr(sce_pseudo, group = sce_pseudo$Subregion)

summary(rowData(sce_pseudo)$high_expr_group_cluster)
#   Mode   FALSE    TRUE 
#logical     426   13416 

## Now filter
sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$high_expr_group_cluster, ]
> dim(sce_pseudo)
# [1] 13416     3

normd <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo))
rownames(normd)<-rownames(sce_pseudo)
colnames(normd)<-sce_pseudo$Subregion

write.table(normd, here::here("processed-data", "09_psychiatric_GSEA","LDSC", "specificity_score_NHP", "NHP_sce_aggregated.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
