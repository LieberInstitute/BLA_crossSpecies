library(dplyr)
library(here)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(tidyr)
library(scater)


# directories
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "08_species_comparisons", "03_DE_analysis")

sce.excit <- readRDS(here("processed-data", "sce_excit_final_subclusters_annotated.rds"))


plotReducedDim(sce.excit, dimred="UMAP", colour_by="DV_axis")


# ======== Pseudobulk analysis ========
# for DE analysis, it makes the most sense to subset to anatomical samples
# and perform pseudobulk and DE analysis. This avoids the statistical "double dipping"
# of using markers genes based on clusters... which are based on DEGs.

# subset SCE to just Later and Basal
sce.subset <- sce.excit[, sce.excit$Subregion %in% c("Lateral", "Basal")]
dim(sce.subset)
# [1] 13874 53306

# aggregate across cells
pseudo <- aggregateAcrossCells(sce.subset, ids = colData(sce.subset)[,c("Subregion", "species", "Sample")])
dim(pseudo)
#[1] 13874    32

unique(pseudo$species)
# [1] "baboon"  "macaque"

# filter low expressing genes
keep <- edgeR::filterByExpr(pseudo, group=pseudo$Subregion)
pseudo <- pseudo[keep,]
dim(pseudo)


# ========= DESeq2 analysis =========

library(DESeq2)

# drop sizeFactors from colData
colData(pseudo)$sizeFactor <- NULL
colnames(colData(pseudo))

dds <- DESeqDataSetFromMatrix(countData=counts(pseudo), 
                              colData=colData(pseudo), 
                              design=~species + Subregion)
dds
# class: DESeqDataSet 
# dim: 10379 86 
# metadata(1): version
# assays(1): counts
# rownames(10379): AURKAIP1 DVL1 ... EMD R3HDM4
# rowData names(0):
#     colnames: NULL
# colData names(66): orig.ident nCount_originalexp ... Sample ncells


dds <- DESeq(dds)

rlogcounts <- rlog(dds)
sampleinfo <- as.data.frame(colData(dds))

# run PCA
library(ggfortify)
pcDat <- prcomp(t(assay(rlogcounts)))

png(here(plot_dir, "Pseudobulk_PCA_PC1vPC2.png"), width=4.5, height=3.5, units="in", res=300)
autoplot(pcDat,
         data = sampleinfo, 
         colour="Subregion", 
         shape="Species",
         size=5,
         x=1,
         y=2) +
    ggtitle("PCA of NHP pseudobulked samples") 
dev.off()


res <- results(dds, contrast=c("Subregion", "Lateral", "Basal"))
res <- lfcShrink(dds,
                 contrast = c("Subregion", "Lateral", "Basal"), res=res, type = 'normal')

res <- as.data.frame(res)
res$gene_name <- rownames(res)

# get log3foldChange > 1 and baseMean > 3000
lateral_top_5 <- res  |> 
    as.data.frame() |> 
    filter(log2FoldChange > 1 & baseMean > 3000) |> 
    arrange(desc(log2FoldChange)) |> 
    head(n=5)
lateral_top_5

# get log3foldChange > 1 and baseMean > 3000
basal_top_5 <- res  |> 
    as.data.frame() |> 
    filter(log2FoldChange < -1 & baseMean > 3000) |> 
    arrange((log2FoldChange)) |> 
    head(n=5)
basal_top_5

# generate a better volcano plot using ggplot
res$sig <- ifelse(res$padj < 0.001 & res$log2FoldChange > 1, "up",
                  ifelse(res$padj < 0.001 & res$log2FoldChange < -1, "down", NA))

png(here(plot_dir, "Volcano_LAvsBA_pseudobulk_ggplot.png"), width=4, height=4, units="in", res=300)
ggplot(res, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color=sig), alpha=0.6) +
    theme_linedraw() +
    theme(legend.position="none") +
    labs(title="LA vs BA",
         #subtitle="ITC pseudobulk DGE",
         x="log2FoldChange",
         y="-log10(p-adj)") +
    geom_text_repel(data=lateral_top_5, aes(x=log2FoldChange, y=-log10(padj), label=gene_name), vjust=1, hjust=1) +
    geom_text_repel(data=basal_top_5, aes(x=log2FoldChange, y=-log10(padj), label=gene_name), vjust=1, hjust=1) +
    xlim(c(-5, 5)) +
    # increase font sizes
    theme(axis.title=element_text(size=16),
          axis.text=element_text(size=12),
          plot.title=element_text(size=16),
          plot.subtitle=element_text(size=14))
dev.off()


# get DEGs < .001 padj and > 1 log2FoldChange from res
res <- res[order(res$padj),]
degs <- res[res$padj < 0.001,]
lateral_degs <- degs[res$log2FoldChange > 1,]
basal_degs <- degs[res$log2FoldChange < -1,]
dim(lateral_degs)
# [1] 146   8

dim(basal_degs)
# [1] 218   8

# save
write.csv(lateral_degs, here("processed-data", "08_species_comparisons", "pseudobulk_lateral_degs.csv"))
write.csv(basal_degs, here("processed-data", "08_species_comparisons", "pseudobulk_basal_degs.csv"))


# ========== Boxplots of average marker gene expression across samples / species ============

# combine putative LA, BA, and aBA clusters
LA_clusters <- c("GULP1_TRHDE", "ZBTB20_SLC4A4")
BA_clusters <- c("PEX5L_MYRIP", "MEIS2_COL25A1")
#aBA_clusters <- c("ESR1_ADRA1A", "GRIK3_TNS3")

sce.excit$subregion_celltype <- NA
sce.excit$subregion_celltype[sce.excit$fine_celltype %in% LA_clusters] <- "LA"
sce.excit$subregion_celltype[sce.excit$fine_celltype %in% BA_clusters] <- "BA"
#sce.excit$subregion_celltype[sce.excit$fine_celltype %in% aBA_clusters] <- "aBA"

# drop NA
sce.subset <- sce.excit[,!is.na(sce.excit$subregion_celltype)]

# aggregate across cells
pseudo <- aggregateAcrossCells(sce.subset, ids = colData(sce.subset)[,c("subregion_celltype", "species", "Sample")])
dim(pseudo)
# [1] 13842 94 

unique(pseudo$species)
# [1] "baboon"  "human"   "macaque"


# drop sizeFactors from colData
colData(pseudo)$sizeFactor <- NULL
colnames(colData(pseudo))

dds <- DESeqDataSetFromMatrix(countData=counts(pseudo), 
                              colData=colData(pseudo), 
                              design=~species + subregion_celltype)
dds

gene2plot <- c("GULP1", "COL25A1")
sampleinfo <- as.data.frame(colData(dds))

rlogcounts <- rlog(dds)
df <- data.frame(assay(rlogcounts))
dim(df)
# [1] 5287   64

df <- df[gene2plot,]
df <- t(df)
df <- cbind(df, sampleinfo)
df <- df[,c(gene2plot, "species", "subregion_celltype","Sample")]

df<- df |> 
    as.data.frame() |> 
    tidyr::pivot_longer(cols=gene2plot, names_to="gene", values_to="expression")

df$gene <- factor(df$gene, levels=gene2plot)

png(here(plot_dir, "Boxplot_putative_LAvBA_marker_genes.png"), width=6, height=4, units="in", res=300)
ggplot(df, aes(x=subregion_celltype, y=expression, fill=subregion_celltype)) +
    geom_boxplot() +
    geom_point(aes(group=subregion_celltype)) +
    geom_line(aes(group = interaction(Sample, gene)), alpha = 0.6) +
    ggh4x::facet_grid2(vars(gene), vars(species)) +
    labs(title="Putative marker genes across species",
         y="Logcounts",
         x=NULL) +
    theme(axis.title=element_text(size=16),
          axis.text=element_text(size=12),
          plot.title=element_text(size=16),
          plot.subtitle=element_text(size=14)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

dev.off()



# =========== Extending this to human samples ============

# "Lateral" clusters
# GULP1_NTNG1
# GULP1_ZBTB20

# "Basal" clusters
# COL25A1_PEX5L
# CHRM3_TRPS1

# "Accessory Basal" clusters
# GPC5_ESR1

# name the fine_type clusters in sce.excit
sce.excit$Subregion_DE <- c()

# Lateral clusters
sce.excit$Subregion_DE[sce.excit$fine_type == "GULP1_NTNG1"] <- "Lateral"
sce.excit$Subregion_DE[sce.excit$fine_type == "GULP1_ZBTB20"] <- "Lateral"

# Basal clusters
sce.excit$Subregion_DE[sce.excit$fine_type == "COL25A1_PEX5L"] <- "Basal"
sce.excit$Subregion_DE[sce.excit$fine_type == "CHRM3_TRPS1"] <- "Basal"

# Accessory Basal clusters
#sce.excit$Subregion_DE[sce.excit$fine_type == "GPC5_ESR1"] <- "Accessory_Basal"

pdf(here(plot_dir, "UMAPs_excitatory_Subregion_DE.pdf"), width = 8, height = 6)
plotReducedDim(sce.excit, dimred = "UMAP", colour_by = "Subregion_DE")
dev.off()


# === heatmap showing DEGS by subregion and species ===

plotGroupedHeatmap(sce.excit, features= features, group="Subregion_DE", block="species", center=TRUE)


library(scuttle)
library(ComplexHeatmap)

# get mean vals
summary <- summarizeAssayByGroup(logcounts(sce.excit),ids = colData(sce.excit)[,c("Subregion_DE", "species")], subset.row=features)
mean.vals <- assays(summary)$mean

# make column annotations: species and Subregion
col_ann <- HeatmapAnnotation(
    subregion=colData(summary)$Subregion_DE,
    species=colData(summary)$species,
    col=list(subregion=c("Basal"="red", "Lateral"="blue"),
             species=c(human = "#AE93BEFF", baboon = "#B4DAE5FF", macaque = "#F0D77BFF"))
)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(7,"RdBu"))(100))

# scale
heat.vals <- t(scale(t(assays(summary)$mean), center=TRUE, scale=TRUE))
pdf(here(plot_dir, "Heatmap_Subregion_DE_top_5_across_species.pdf"), width = 6, height = 5)
ComplexHeatmap::Heatmap(heat.vals,
                        top_annotation=col_ann,
                        col=cols,
                        name="mean logcounts",
                        )
dev.off()