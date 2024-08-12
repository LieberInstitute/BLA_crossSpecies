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

lateral_degs <- read.csv(here("processed-data", "08_species_comparisons", "pseudobulk_lateral_degs.csv"))
basal_degs <- read.csv(here("processed-data", "08_species_comparisons", "pseudobulk_basal_degs.csv"))



# ===== compare clusters =====

gene_list <- list(Lateral = lateral_degs$gene_name, Basal = basal_degs$gene_name)

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


simplified_cc <- simplify(ego_compare_cc, cutoff=0.7, by="p.adjust", select_fun=min)
simplified_bp <- simplify(ego_compare_bp, cutoff=0.7, by="p.adjust", select_fun=min)
simplified_mf <- simplify(ego_compare_mf, cutoff=0.7, by="p.adjust", select_fun=min)


png(here(plot_dir, "LAvBA_CC_term_simplified.png"), width=4.75, height=4.75, units="in", res=300)
dotplot(simplified_cc, showCategory = 5, title = "Cellular Components") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x="Subregion")
dev.off()

png(here(plot_dir, "LAvBA_BP_term_simplified.png"), width=5, height=8, units="in", res=300)
dotplot(simplified_bp, showCategory = 5, title = "Biological Process") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x="Subregion")
dev.off()

png(here(plot_dir, "LAvBA_MF_term_simplified.png"), width=5, height=8, units="in", res=300)
dotplot(simplified_mf, showCategory = 5, title = "Molecular Function")  +
    labs(x="Subregion")
dev.off()



# ======= Get genes from GO terms ========

# load excitatory SCE data and subset to putative LA vs BA clusters
sce.excit <- readRDS(here("processed-data", "sce_excit_final_subclusters_annotated.rds"))
sce.subset <- sce.excit[, sce.excit$fine_celltype %in% c("GULP1_TRHDE", "ZBTB20_SLC4A4", "PEX5L_MYRIP", "MEIS2_COL25A1")]
sce.subset

# Get BP terms
specific_terms <- simplified_bp@compareClusterResult  %>%
    filter(Description %in% c("response to nicotine", "cell-cell adhesion via plasma-membrane adhesion molecules"))

nicotine_genes <- unlist(strsplit(specific_terms$geneID[specific_terms$Description == "response to nicotine"], "/"))
adhesion_genes <- unlist(strsplit(specific_terms$geneID[specific_terms$Description == "cell-cell adhesion via plasma-membrane adhesion molecules"], "/"))
bp_genes <- c(nicotine_genes[1:5], adhesion_genes[1:5])

# Get CC Terms
specific_terms <- simplified_cc@compareClusterResult  %>%
    filter(Description %in% c("synaptic membrane", "cation channel complex"))

synaptic_genes <- unlist(strsplit(specific_terms$geneID[specific_terms$Description == "synaptic membrane"], "/"))
channel_genes <- unlist(strsplit(specific_terms$geneID[specific_terms$Description == "cation channel complex"], "/"))
cc_genes <- c(synaptic_genes[1:5], channel_genes[1:5])

# get MF terms
specific_terms <- simplified_mf@compareClusterResult  %>%
    filter(Description %in% c("neurotransmitter receptor activity", "calcium channel activity"))

receptor_genes <- unlist(strsplit(specific_terms$geneID[specific_terms$Description == "neurotransmitter receptor activity"], "/"))
calcium_genes <- unlist(strsplit(specific_terms$geneID[specific_terms$Description == "calcium channel activity"], "/"))
mf_genes <- c(receptor_genes[1:5], calcium_genes[1:5])


# === Grouped heatmap of genes ===
plotGroupedHeatmap(sce.subset, features=bp_genes, group="fine_celltype", block="species",center=TRUE, scale=TRUE)
plotGroupedHeatmap(sce.subset, features=cc_genes, group="fine_celltype", block="species",center=TRUE, scale=TRUE)
plotGroupedHeatmap(sce.subset, features=mf_genes, group="fine_celltype", block="species",center=TRUE, scale=TRUE)

