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

sce.inhib <- readRDS(here(processed_dir, "sce.inhib.final.rds"))
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
                              design=~DV_axis)
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
         x=3,
         y=4) +
    ggtitle("PCs 6 seperates DV axis") 
dev.off()


plotPCA(rlogcounts, intgroup=c("DV_axis"), ntop=2000)


# DE analysis

res <- results(dds, contrast=c("Subregion", "Lateral", "Basal"))
res <- lfcShrink(dds,
                 contrast = c("Subregion", "Lateral", "Basal"), res=res, type = 'normal')


EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                #FCcutoff = 1,
                pCutoff = 10e-6,
                #selectLab=c('COL25A1'),
                title = 'Lateral vs Basal',
                subtitle='pseudobulk analysis'
)



# ======= DV axis DE analysis =======
# subset to macaque, then lateral and basal
sce.subset <- sce.excit[, sce.excit$species == "macaque"]
sce.subset <- sce.subset[, sce.subset$DV_axis %in% c("Dorsal", "Ventral")]
sce.subset$DV_axis <- as.factor(sce.subset$DV_axis)

sce.lateral <- sce.subset[, sce.subset$Subregion == "Lateral"]
sce.basal <- sce.subset[, sce.subset$Subregion == "Basal"]

# pseudobulk
pb.lateral <- aggregateAcrossCells(sce.lateral, ids = colData(sce.lateral)[,c("DV_axis", "Sample")])
dim(pb.lateral)

pb.basal <- aggregateAcrossCells(sce.basal, ids = colData(sce.basal)[,c("DV_axis", "Sample")])
dim(pb.basal)

# filter low expressing genes
keep <- edgeR::filterByExpr(pb.lateral, group=pb.lateral$DV_axis)
pb.lateral <- pb.lateral[keep,]
dim(pb.lateral)

keep <- edgeR::filterByExpr(pb.basal, group=pb.basal$DV_axis)
pb.basal <- pb.basal[keep,]
dim(pb.basal)

# ===== DESeq objects =====
# lateral
colData(pb.lateral)$sizeFactor <- NULL
dds.lateral <- DESeqDataSetFromMatrix(countData=counts(pb.lateral), 
                              colData=colData(pb.lateral), 
                              design=~DV_axis)
dds.lateral

dds.lateral <- DESeq(dds.lateral)


rlogcounts <- rlog(dds.lateral)
sampleinfo <- as.data.frame(colData(dds.lateral))

# DE analysis
res.lateral <- results(dds.lateral, contrast=c("DV_axis", "Dorsal", "Ventral"))

EnhancedVolcano(res.lateral,
                lab = rownames(res.lateral),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1.5,
                pCutoff = .05,
                #selectLab=c('PEX5L'),
                title = 'Lateral: Dorsal vs Ventral',
                subtitle='PB analysis of macaque samples'
)


# basal
colData(pb.basal)$sizeFactor <- NULL
dds.basal <- DESeqDataSetFromMatrix(countData=counts(pb.basal), 
                                      colData=colData(pb.basal), 
                                      design=~DV_axis)
dds.basal

dds.basal <- DESeq(dds.basal)


rlogcounts <- rlog(dds.basal)
sampleinfo <- as.data.frame(colData(dds.basal))

# DE analysis
res.basal <- results(dds.basal, contrast=c("DV_axis", "Dorsal", "Ventral"))

EnhancedVolcano(res.basal,
                lab = rownames(res.basal),
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1.5,
                pCutoff = .01,
                #selectLab=c('PEX5L'),
                title = 'Basal: Dorsal vs Ventral',
                subtitle='PB analysis of macaque samples'
)
 