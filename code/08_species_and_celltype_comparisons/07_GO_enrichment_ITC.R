library(here)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(tidyr)
library(scater)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)

# directories
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "08_species_comparisons")

itc1_degs <- read.csv(here("processed-data", "08_species_comparisons", "pseudobulk_itc1_degs.csv"))
itc2_degs <- read.csv(here("processed-data", "08_species_comparisons", "pseudobulk_itc2_degs.csv"))



# ===== compare clusters =====

gene_list <- list(TSHZ1.1 = itc1_degs$gene_name, TSHZ1.2 = itc2_degs$gene_name)

# Perform GO enrichment analysis using compareCluster
ego_compare_cc <- compareCluster(geneCluster = gene_list, 
                              fun = "enrichGO", 
                              OrgDb = org.Hs.eg.db,
                              keyType = "SYMBOL",
                              ont = "CC", 
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.05,
                              readable = TRUE)

ego_compare_bp <- compareCluster(geneCluster = gene_list,
                              fun = "enrichGO",
                              OrgDb = org.Hs.eg.db,
                              keyType = "SYMBOL",
                              ont = "BP",
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.05,
                              readable = TRUE)

ego_compare_mf <- compareCluster(geneCluster = gene_list,
                                 fun = "enrichGO",
                                 OrgDb = org.Hs.eg.db,
                                 keyType = "SYMBOL",
                                 ont = "MF",
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.05,
                                 readable = TRUE)


png(here(plot_dir, "itc_GO_term_comparisona.png"), width=14, height=12, units="in", res=300)

p1 <- dotplot(ego_compare_cc, showCategory = 20, title = "Cellular Components") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- dotplot(ego_compare_bp, showCategory = 10, title = "Biological Process") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- dotplot(ego_compare_mf, showCategory = 20, title = "Molecular Function") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1+p2+p3
dev.off()
