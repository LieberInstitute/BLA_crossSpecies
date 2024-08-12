library(dplyr)
library(here)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(scater)
library(spatialLIBD)

# directories
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "08_species_comparisons")

sce.inhib <- readRDS(here("processed-data", "sce_inhib_final_subclusters_annotated.rds"))
sce.inhib

rowData(sce.inhib)$gene_name <- rownames(sce.inhib)


# ======== PB ITC clusters ========

# drop all other cell-types not containing "TSHZ1"
sce.itc <- sce.inhib[,grepl("TSHZ1", colData(sce.inhib)$fine_celltype)]
sce.itc

pb.itc <- aggregateAcrossCells(sce.itc, 
                               id=colData(sce.itc)[,c("fine_celltype", "Sample")])
pb.itc

# factorize fine_celltypes
pb.itc$fine_celltype <- factor(pb.itc$fine_celltype)

# Removing all pb.itc-bulk samples with 'insufficient' cells.
pb.itc <- pb.itc[,pb.itc$ncells >= 10]
pb.itc

keep <- edgeR::filterByExpr(pb.itc, group=pb.itc$fine_celltype)
pb.itc <- pb.itc[keep,]
dim(pb.itc)
# [1] 5287   64


# ========= DESeq2 analysis =========

library(DESeq2)

# drop sizeFactors from colData
colData(pb.itc)$sizeFactor <- NULL
colnames(colData(pb.itc))

dds <- DESeqDataSetFromMatrix(countData=counts(pb.itc), 
                              colData=colData(pb.itc), 
                              design=~species + fine_celltype)
dds
# class: DESeqDataSet 
# dim: 5287 64 
# metadata(1): version
# assays(1): counts
# rownames(5287): ACAP3 NADK ... CUL2 ABLIM2
# rowData names(0):
#     colnames: NULL
# colData names(39): orig.ident nCount_originalexp ... ncells itc


# === PCA ===
library(ggfortify)
rlogcounts <- rlog(dds)
sampleinfo <- as.data.frame(colData(dds))

pcDat <- prcomp(t(assay(rlogcounts)))

png(here(plot_dir, "Pseudobulk_PCA_PC1vPC2_ITCs.png"), width=6.5, height=5, units="in", res=300)
autoplot(pcDat,
         data = sampleinfo, 
         colour="species", 
         shape="fine_celltype",
         size=5,
         x=1,
         y=2) +
    ggtitle("PCs of pseudobulk ITCs across species") 
dev.off()

# === DESeq2 ===

dds <- DESeq(dds)

res <- results(dds, contrast=c("fine_celltype", "TSHZ1.1", "TSHZ1.2"))
res <- lfcShrink(dds,
                 contrast = c("fine_celltype", "TSHZ1.1", "TSHZ1.2"), res=res, type = 'normal')
res <- as.data.frame(res)
res$gene_name <- rownames(res)

# get log2foldChange > 1 and baseMean > 25
itc1_top_5 <- res  |> 
    as.data.frame() |> 
    filter(log2FoldChange > 1 & baseMean > 25) |> 
    arrange(desc(log2FoldChange)) |> 
    head(n=5)
itc1_top_5

# get log2foldChange > 1 and baseMean > 25
itc2_top_5 <- res  |> 
    as.data.frame() |> 
    filter(log2FoldChange < -1 & baseMean > 25) |> 
    arrange((log2FoldChange)) |> 
    head(n=5)
itc2_top_5

# generate a better volcano plot using ggplot
res$sig <- ifelse(res$padj < 0.0001 & res$log2FoldChange > 1, "up",
                  ifelse(res$padj < 0.0001 & res$log2FoldChange < -1, "down", NA))

png(here(plot_dir, "Volcano_ITCs_psuedobulk_ggplot.png"), width=5, height=5, units="in", res=300)
ggplot(res, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color=sig), alpha=0.6) +
    theme_linedraw() +
    theme(legend.position="none") +
    labs(title="TSHZ1.1 vs TSHZ1.2",
         #subtitle="ITC pseudobulk DGE",
         x="log2FoldChange",
         y="-log10(p-adj)") +
    geom_text_repel(data=itc1_top_5, aes(x=log2FoldChange, y=-log10(padj), label=gene_name), vjust=1, hjust=1) +
    geom_text_repel(data=itc2_top_5, aes(x=log2FoldChange, y=-log10(padj), label=gene_name), vjust=1, hjust=1) +
    xlim(c(-5, 5)) +
    # increase font sizes
    theme(axis.title=element_text(size=16),
          axis.text=element_text(size=12),
          plot.title=element_text(size=16),
          plot.subtitle=element_text(size=14))
dev.off()


# get DEGs < .001 padj and > 1 log2FoldChange from res
res <- res[order(res$padj),]
degs <- res[res$padj < 0.0001,]
itc1_degs <- degs[res$log2FoldChange > 1,]
itc2_degs <- degs[res$log2FoldChange < -1,]
dim(itc1_degs)
#[1] 105   8

dim(itc2_degs)
# [1] 116   8

# save
write.csv(itc1_degs, here("processed-data", "08_species_comparisons", "pseudobulk_itc1_degs.csv"))
write.csv(itc2_degs, here("processed-data", "08_species_comparisons", "pseudobulk_itc2_degs.csv"))



# ========= Boxplots average expression per sample, split by species ========

top5_genes <- c(itc1_top_5$gene_name, itc2_top_5$gene_name)
sampleinfo <- as.data.frame(colData(dds))

df <- data.frame(assay(rlogcounts))
dim(df)
# [1] 5287   64

df <- df[top5_genes,]
df <- t(df)
df <- cbind(df, sampleinfo)
df <- df[,c(top5_genes, "species", "fine_celltype","Sample")]

df_itc1 <- df |> 
    as.data.frame() |> 
    tidyr::pivot_longer(cols=1:5, names_to="gene", values_to="expression")

df_itc2 <- df |> 
    as.data.frame() |> 
    tidyr::pivot_longer(cols=6:10, names_to="gene", values_to="expression")

df_itc1$gene <- factor(df_itc1$gene, levels=top5_genes)
df_itc2$gene <- factor(df_itc2$gene, levels=top5_genes)

png(here(plot_dir, "Boxplot_top5_ITC1_genes.png"), width=7, height=5, units="in", res=300)
ggplot(df_itc1, aes(x=fine_celltype, y=expression, fill=fine_celltype)) +
    geom_boxplot() +
    geom_point(aes(group=fine_celltype)) +
    geom_line(aes(group = interaction(Sample, gene)), alpha = 0.6) +
    ggh4x::facet_grid2(vars(species), vars(gene)) +
    labs(title="Top 5 DE genes in ITC1",
         x="Gene",
         y="Expression") +
    theme(axis.title=element_text(size=16),
          axis.text=element_text(size=12),
          plot.title=element_text(size=16),
          plot.subtitle=element_text(size=14)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
dev.off()


png(here(plot_dir, "Boxplot_top5_ITC2_genes.png"), width=7, height=5, units="in", res=300)
ggplot(df_itc2, aes(x=fine_celltype, y=expression, fill=fine_celltype)) +
    geom_boxplot() +
    geom_point(aes(group=fine_celltype)) +
    geom_line(aes(group = interaction(Sample, gene)), alpha = 0.6) +
    ggh4x::facet_grid2(vars(species), vars(gene)) +
    labs(title="Top 5 DE genes in ITC2",
         x="Gene",
         y="Expression") +
    theme(axis.title=element_text(size=16),
          axis.text=element_text(size=12),
          plot.title=element_text(size=16),
          plot.subtitle=element_text(size=14)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
dev.off()


# ==== heatmap of top 5 genes in ITC1 and ITC2

library(scuttle)
library(ComplexHeatmap)

# get mean vals
features <- top5_genes
summary <- summarizeAssayByGroup(assay(rlogcounts),ids = colData(rlogcounts)[,c("fine_celltype", "species")], subset.row=features)
mean.vals <- assays(summary)$mean

# make column annotations: species and Subregion
celltype_colors <- c("TSHZ1.1"="red", "TSHZ1.2"="blue")
species_colors <- c(human = "#AE93BEFF", baboon = "#B4DAE5FF", macaque = "#F0D77BFF")

col_ann <- HeatmapAnnotation(
    fine_celltype=colData(summary)$fine_celltype,
    species=colData(summary)$species,
    col=list(fine_celltype=c("TSHZ1.1"="red", "TSHZ1.2"="blue"),
             species=species_colors)
)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(7,"RdBu"))(100))

# scale
heat.vals <- t(scale(t(assays(summary)$mean), center=TRUE, scale=TRUE))
png(here(plot_dir, "Heatmap_Subregion_DE_top_5_across_species.png"), width = 5, height = 3, units="in", res=300)
ComplexHeatmap::Heatmap(heat.vals,
                        top_annotation=col_ann,
                        col=cols,
                        name="centered.scaled"
)


dev.off()

