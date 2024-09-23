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
processed_dir = here("processed-data", "07_annotation_and_characterization")
plot_dir = here("plots", "08_species_comparisons", "05_DE_PRKCD")

sce.inhib <- readRDS(here("processed-data", "07_annotation_and_characterization", "sce_inhib_final_subclusters_annotated.rds"))
sce.inhib

rowData(sce.inhib)$gene_name <- rownames(sce.inhib)


# ======== PB prkcd clusters ========

# drop all other cell-types not containing "TSHZ1"
sce.prkcd <- sce.inhib[,grepl("PRKCD", colData(sce.inhib)$fine_celltype)]
sce.prkcd

pb.prkcd <- aggregateAcrossCells(sce.prkcd, 
                               id=colData(sce.prkcd)[,c("fine_celltype", "Sample")])
pb.prkcd

# factorize fine_celltypes
pb.prkcd$fine_celltype <- factor(pb.prkcd$fine_celltype)

# Removing all pb.prkcd-bulk samples with 'insufficient' cells.
pb.prkcd <- pb.prkcd[,pb.prkcd$ncells >= 10]
pb.prkcd

keep <- edgeR::filterByExpr(pb.prkcd, group=pb.prkcd$fine_celltype)
pb.prkcd <- pb.prkcd[keep,]
dim(pb.prkcd)
# [1] 7118   18

# drop human samples
pb.prkcd <- pb.prkcd[,pb.prkcd$species == c("macaque")]
dim(pb.prkcd)
# [1] 7118   14

# ========= DESeq2 analysis =========

library(DESeq2)

# drop sizeFactors from colData
colData(pb.prkcd)$sizeFactor <- NULL
colnames(colData(pb.prkcd))

dds <- DESeqDataSetFromMatrix(countData=counts(pb.prkcd), 
                              colData=colData(pb.prkcd), 
                              design=~fine_celltype)
dds
# class: DESeqDataSet 
# dim: 5287 64 
# metadata(1): version
# assays(1): counts
# rownames(5287): ACAP3 NADK ... CUL2 ABLIM2
# rowData names(0):
#     colnames: NULL
# colData names(39): orig.ident nCount_originalexp ... ncells prkcd


# === PCA ===
library(ggfortify)
rlogcounts <- rlog(dds)
sampleinfo <- as.data.frame(colData(dds))

pcDat <- prcomp(t(assay(rlogcounts)))

png(here(plot_dir, "Pseudobulk_PCA_PC1vPC2_prkcds.png"), width=6.5, height=5, units="in", res=300)
autoplot(pcDat,
         data = sampleinfo, 
         colour="species", 
         shape="fine_celltype",
         size=5,
         x=1,
         y=2) +
    ggtitle("PCs of pseudobulk prkcds across species") 
dev.off()

# === DESeq2 ===

dds <- DESeq(dds)

res <- results(dds, contrast=c("fine_celltype", "PRKCD_DRD1", "PRKCD_DRD2"))
res <- lfcShrink(dds,
                 contrast = c("fine_celltype",  "PRKCD_DRD1", "PRKCD_DRD2"), res=res, type = 'normal')
res <- as.data.frame(res)
res$gene_name <- rownames(res)

# get log2foldChange > 1 and baseMean > 25
prkcd1_top_5 <- res  |> 
    as.data.frame() |> 
    filter(log2FoldChange > 1) |> 
    arrange(desc(log2FoldChange)) |> 
    head(n=5)
prkcd1_top_5

# get log2foldChange > 1 and baseMean > 25
prkcd2_top_5 <- res  |> 
    as.data.frame() |> 
    filter(log2FoldChange < -1) |> 
    arrange((log2FoldChange)) |> 
    head(n=5)
prkcd2_top_5

# generate a better volcano plot using ggplot
res$sig <- ifelse(res$padj < 0.05 & res$log2FoldChange > 1, "up",
                  ifelse(res$padj < 0.05 & res$log2FoldChange < -1, "down", NA))

png(here(plot_dir, "Volcano_prkcds_psuedobulk_ggplot.png"), width=5, height=5, units="in", res=300)
ggplot(res, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color=sig), alpha=0.6) +
    theme_linedraw() +
    theme(legend.position="none") +
    labs(title= "PRKCD_DRD1", "PRKCD_DRD2",
         #subtitle="prkcd pseudobulk DGE",
         x="log2FoldChange",
         y="-log10(p-adj)") +
    geom_text_repel(data=subset(res, gene_name %in% c("DRD3", "GRIK2", "PDYN", "EBF1", "GRIK4", "CNR1", "FOXP2", "ERBB4", "TAC1", "SST")),
                aes(x=log2FoldChange, y=-log10(padj), label=gene_name),  vjust=1, hjust=-.75) +
    geom_text_repel(data=subset(res, gene_name %in% c("DRD2", "GRIK3",  "PENK", "EYA4", "CHRM3", "PTPRM", "HTR2C", "CALCRL","HCN1")),
                aes(x=log2FoldChange, y=-log10(padj), label=gene_name),  vjust=1, hjust=1.25) +
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
prkcd1_degs <- degs[res$log2FoldChange > 1,]
prkcd2_degs <- degs[res$log2FoldChange < -1,]
dim(prkcd1_degs)
#[1] 105   8

dim(prkcd2_degs)
# [1] 116   8

# save
write.csv(prkcd1_degs, here("processed-data", "08_species_comparisons","05_DE_prkcds", "pseudobulk_prkcd_drd1_degs.csv"))
write.csv(prkcd2_degs, here("processed-data", "08_species_comparisons","05_DE_prkcds", "pseudobulk_prkcd_drd2_degs.csv"))



# ========= Boxplots average expression per sample, split by species ========

prkcd1_genes2plot <- c("DRD3", "GRIK2", "PDYN", "GRIK4", "CNR1", "TAC1", "SST")
prkcd2_genes2plot <- c("DRD2", "GRIK3",  "PENK",  "CHRM3", "HTR2C", "CALCRL")

top5_genes <- c(prkcd1_genes2plot, prkcd2_genes2plot)
sampleinfo <- as.data.frame(colData(dds))

df <- data.frame(assay(rlogcounts))
dim(df)
# [1] 5287   64

df <- df[top5_genes,]
df <- t(df)
df <- cbind(df, sampleinfo)
df <- df[,c(top5_genes, "fine_celltype","Sample")]

df_prkcd1 <- df |> 
    as.data.frame() |> 
    tidyr::pivot_longer(cols=1:7, names_to="gene", values_to="expression")

df_prkcd2 <- df |> 
    as.data.frame() |> 
    tidyr::pivot_longer(cols=8:13, names_to="gene", values_to="expression")

df_prkcd1$gene <- factor(df_prkcd1$gene, levels=top5_genes)
df_prkcd2$gene <- factor(df_prkcd2$gene, levels=top5_genes)

png(here(plot_dir, "Boxplot_top5_prkcd_drd1_genes.png"), width=7, height=5, units="in", res=300)
ggplot(df_prkcd1, aes(x=fine_celltype, y=expression, fill=fine_celltype)) +
    geom_boxplot() +
    geom_point(aes(group=fine_celltype)) +
    geom_line(aes(group = interaction(Sample, gene)), alpha = 0.6) +
    facet_grid(cols=vars(gene)) +
    labs(title="Pseudobulk DE genes in PRKCD_DRD1",
         x="Gene",
         y="Expression") +
    theme(axis.title=element_text(size=16),
          axis.text=element_text(size=12),
          plot.title=element_text(size=16),
          plot.subtitle=element_text(size=14)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
dev.off()


png(here(plot_dir, "Boxplot_top5_prkcd_drd2_genes.png"), width=7, height=5, units="in", res=300)
ggplot(df_prkcd2, aes(x=fine_celltype, y=expression, fill=fine_celltype)) +
    geom_boxplot() +
    geom_point(aes(group=fine_celltype)) +
    geom_line(aes(group = interaction(Sample, gene)), alpha = 0.6) +
    facet_grid(cols=vars(gene)) +
    labs(title="Pseudobulk DE genes in PRKCD_DRD2",
         x="Gene",
         y="Expression") +
    theme(axis.title=element_text(size=16),
          axis.text=element_text(size=12),
          plot.title=element_text(size=16),
          plot.subtitle=element_text(size=14)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
dev.off()
