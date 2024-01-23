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

res <- results(dds, contrast=c("Subregion", "Lateral", "Basal"))
res <- lfcShrink(dds,
                 contrast = c("Subregion", "Lateral", "Basal"), res=res, type = 'normal')

pdf(here(plot_dir, "Volcano_Lateral_vs_Basal_pseudo.pdf"), width = 7, height = 8)
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
dev.off()



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


features <- c(rownames(lateral_top_5), rownames(basal_top_5))
pdf(here(plot_dir, "Violins_Subregion_DE_top_5_pseudo.pdf"), width = 6, height = 2.5)
plotExpression(sce.subset, features = features, colour_by = "Subregion", x="Subregion", ncol=5) +
    theme(legend.position="none")
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

# === UMAPs with top marker genes ===

library(escheR)

rowData(sce.excit)$gene_name <- rownames(sce.excit)
sce.temp <- sce.excit
sce.temp$expr_GULP1 <- logcounts(sce.excit)[which(rowData(sce.excit)$gene_name=="GULP1"),]
sce.temp$expr_COL25A1 <- logcounts(sce.excit)[which(rowData(sce.excit)$gene_name=="COL25A1"),]

# Point Binning version
p <- make_escheR(sce.temp, dimred='UMAP') 

p1 <- p |>
        add_ground_bin(var="Subregion_DE") |> 
        add_fill_bin(var = "expr_GULP1") +
        # Customize aesthetics
        scale_fill_gradient(low = "white", high = "black", name = "GULP1")+
        theme_void()

p2<- p |>
        add_ground_bin(var="Subregion_DE") |> 
        add_fill_bin(var = "expr_COL25A1") +
        # Customize aesthetics
        scale_fill_gradient(low = "white", high = "black", name = "COL25A1")+
        theme_void()

pdf(here(plot_dir, "UMAPs_GULP1_COL25A1_escheR.pdf"), width = 8, height = 4)
p1+p2
dev.off()

pdf(here(plot_dir, "UMAPs_GULP1_COL25A1.pdf"), width = 8, height = 8)
p1 <- plotReducedDim(sce.temp, dimred = "UMAP", colour_by = "Subregion", point_size=0.1) +
    ggtitle("Subregion Samples") +
    theme_void()
p2 <- plotReducedDim(sce.temp, dimred = "UMAP", colour_by = "fine_type", point_size=0.1) +
    ggtitle("Clusters") +
    theme(legend.position="none") +
    theme_void()
p3 <- plotReducedDim(sce.temp, dimred = "UMAP", colour_by = "GULP1", point_size=0.1) +
    scale_colour_gradient(low = "grey", high = "red", name = "logcount") +
    ggtitle("GULP1 expression") +
    theme_void()
p4 <- plotReducedDim(sce.temp, dimred = "UMAP", colour_by = "COL25A1", point_size=0.1) +
    scale_colour_gradient(low = "grey", high = "red", name = "logcount") +
    ggtitle("COL25A1: expression") +
    theme_void()
(p1+p2)/(p3+p4)
dev.off()



