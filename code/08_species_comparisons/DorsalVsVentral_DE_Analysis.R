library(dplyr)
library(here)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(tidyr)
library(scater)
library(SummarizedExperiment)


# directories
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "08_species_comparisons", "03_DE_analysis")

sce.excit <- readRDS(here(processed_dir, "sce.excit.integrated.annotated.rds"))


# ======== Pseudobulk analysis ========
# for DE analysis, it makes the most sense to subset to anatomical samples
# and perform pseudobulk and DE analysis. This avoids the statistical "double dipping"
# of using markers genes based on clusters... which are based on DEGs.

# subset SCE to just Later and Basal
sce.subset <- sce.excit[, sce.excit$Subregion %in% c("Lateral", "Basal")]
dim(sce.subset)
# [1] 13874 53306

# aggregate across cells
sce.subset$DV_axis <- as.factor(sce.subset$DV_axis)
pseudo <- aggregateAcrossCells(sce.subset, ids = colData(sce.subset)[,c("Subregion", "species","DV_axis", "Sample")])
dim(pseudo)
#[1] 13874    32

unique(pseudo$species)
# [1] "baboon"  "macaque"

# filter low expressing genes
keep <- edgeR::filterByExpr(pseudo, group=pseudo$Subregion)
pseudo <- pseudo[keep,]
dim(pseudo)

library(DESeq2)
# drop sizeFactors from colData
colData(pseudo)$sizeFactor <- NULL
colnames(colData(pseudo))

dds <- DESeqDataSetFromMatrix(countData=counts(pseudo), 
                              colData=colData(pseudo), 
                              design=~species+Subregion)
dds

dds <- DESeq(dds)


rlogcounts <- rlog(dds)
sampleinfo <- as.data.frame(colData(dds))

# run PCA
library(ggfortify)
pcDat <- prcomp(t(assay(rlogcounts)))

png(here(plot_dir, "Pseudobulk_PCA_PC1vPC2.png"), width=6.5, height=5, units="in", res=300)
autoplot(pcDat,
         data = sampleinfo, 
         colour="Subregion", 
         shape="Species",
         size=5,
         x=1,
         y=2) +
    ggtitle("PCs 1 and 2 seperates NHP species and subregions") 
dev.off()

png(here(plot_dir, "Pseudobulk_PCA_PC2vPC6.png"), width=6.5, height=5, units="in", res=300)
autoplot(pcDat,
         data = sampleinfo, 
         colour="DV_axis",
         size=5,
         x=2,
         y=6) +
    ggtitle("PCs 6 seperates DV axis") 
dev.off()

# DE analysis

res <- results(dds, contrast=c("Subregion", "Lateral", "Basal"))
res <- lfcShrink(dds,
                 contrast = c("Subregion", "Lateral", "Basal"), res=res, type = 'normal')


png(here(plot_dir, "Volcano_LateralVsBasal.png"), width=7, height=8, units="in", res=300)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                #FCcutoff = 1,
                pCutoff = 10e-6,
                selectLab=c('GULP1', 'PEX5L',"COL25A1","SATB1"),
                title = 'Lateral vs Basal',
                subtitle='pseudobulk analysis'
)
dev.off()


# get DEGs < .05 padj and > 1 log2FoldChange from res
res <- res[order(res$padj),]
degs <- res[res$padj < 0.05,]
lateral_degs <- degs[res$log2FoldChange > 1,]
basal_degs <- degs[res$log2FoldChange < -1,]
length(lateral_degs)
# 146

length(basal_degs)
# 218

# save
write.csv(lateral_degs, here("processed-data", "08_species_comparisons", "pseudobulk_ateral_degs.csv"))
write.csv(basal_degs, here("processed-data", "08_species_comparisons", "pseudobulk_basal_degs.csv"))

# ======= GO analysis on lateral vs basal DEGs =========

library(clusterProfiler)
library(org.Hs.eg.db)
 
# == Lateral GO enrichment ==
ego.bp <- enrichGO(gene = rownames(lateral_degs),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",  # Use appropriate keyType for your gene list
                ont = "BP",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05)

ego.cc <- enrichGO(gene = rownames(lateral_degs),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",  # Use appropriate keyType for your gene list
                ont = "CC",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05)

ego.mf <- enrichGO(gene = rownames(lateral_degs),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",  # Use appropriate keyType for your gene list
                ont = "MF",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05)


png(here(plot_dir, "Lateral_GO_BiologicalProcess.png"), width=5, height=5, units="in", res=300)
dotplot(ego.bp, showCategory=30) + ggtitle("Lateral GO Biological Process")
dev.off()

png(here(plot_dir, "Lateral_GO_CellularComponent.png"), width=5, height=5, units="in", res=300)
dotplot(ego.cc, showCategory=30) + ggtitle("Lateral GO Cellular Component")
dev.off()

png(here(plot_dir, "Lateral_GO_MolecularFunction.png"), width=5, height=5, units="in", res=300)
dotplot(ego.mf, showCategory=30) + ggtitle("Lateral GO Molecular Function")
dev.off()

# == Basal GO enrichment ==
ego.bp <- enrichGO(gene = rownames(basal_degs),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",  # Use appropriate keyType for your gene list
                ont = "BP",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05)

ego.cc <- enrichGO(gene = rownames(basal_degs),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",  # Use appropriate keyType for your gene list
                ont = "CC",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05)

ego.mf <- enrichGO(gene = rownames(basal_degs),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",  # Use appropriate keyType for your gene list
                ont = "MF",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05)


png(here(plot_dir, "Basal_GO_BiologicalProcess.png"), width=5, height=7.5, units="in", res=300)
dotplot(ego.bp, showCategory=30) + ggtitle("Basal GO Biological Process")
dev.off()

# 0 enriched terms found
#png(here(plot_dir, "Basal_GO_CellularComponent.png"), width=5, height=5, units="in", res=300)
#dotplot(ego.cc, showCategory=30) + ggtitle("Basal GO Cellular Component")
#dev.off()

png(here(plot_dir, "Basal_GO_MolecularFunction.png"), width=5, height=5, units="in", res=300)
dotplot(ego.mf, showCategory=30) + ggtitle("Basal GO Molecular Function")
dev.off()