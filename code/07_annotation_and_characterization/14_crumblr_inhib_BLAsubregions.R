library(SpatialExperiment)
library(sessioninfo)
library(ggplot2)
library(here)
library(crumblr)
library(variancePartition)
library(dreamlet)
library(patchwork)
library(ggbeeswarm)

## save directories
plot_dir = here("plots", "07_annotation_and_characterization","14_crumblr_inhib_BLAsubregions")
processed_dir <- here("processed-data")

# load sce
sce <- readRDS(here("processed-data","05_batch_correction", "within_species", "macaque", "sce_macaque_inhib_label_transfer.rds"))
sce

colnames(colData(sce))
#  [1] "orig.ident"            "nCount_originalexp"    "nFeature_originalexp" 
#  [4] "nCount_SCT"            "nFeature_SCT"          "Sample"               
#  [7] "Barcode"               "sum"                   "detected"             
# [10] "subsets_Mito_sum"      "subsets_Mito_detected" "subsets_Mito_percent" 
# [13] "total"                 "key"                   "subject"              
# [16] "species"               "subregion"             "dv_axis"              
# [19] "sample_num"            "high_mito"             "low_lib"              
# [22] "low_genes"             "discard_auto"          "discard_minimum"      
# [25] "discard"               "doubletScore"          "sex.x"                
# [28] "sizeFactor"            "SCT_snn_res.0.8"       "seurat_clusters"      
# [31] "sex"                   "SCT_snn_res.0.5"       "broad_celltype"       
# [34] "ident"                 "SCT_snn_res.0.7"       "fine_celltype" 


# ========== Dropping CeA ===========

# drop
sce <- sce[, sce$subregion != "Central Nucleus"]

# rename Latera to LA, Basal to BA, and Accessory Basal to aBA. Necessary for later contrasts
sce$subregion <- factor(sce$subregion,
  levels = c("Lateral", "Basal", "Accessory Basal"),
  labels = c("LA", "BA", "aBA")
)


# drop levels
sce$subregion <- as.factor(sce$subregion)
sce$subregion <- droplevels(sce$subregion)

# drop bonafide CeA cell types - "DLK1_ZFHX3", "PENK_DRD2", "TAC1_PPP1R1B", "SST_TAC1"
sce <- sce[,!sce$fine_celltype %in% c("DLK1_ZFHX3", "PENK_DRD2", "TAC1_PPP1R1B", "SST_TAC1")]



# ===========================================================================
#  Creating pseudobulk data for crumblr analysis 
# ===========================================================================

# Create pseudobulk data by specifying cluster_id and sample_id for aggregating cells
pb <- aggregateToPseudoBulk(sce,
  assay = "counts",
  cluster_id = "fine_celltype",
  sample_id = "Sample",
  verbose = FALSE
)

cellCounts(pb)

cobj <- crumblr(cellCounts(pb))

# reorder levels so Lateral is first
pb$subregion <- factor(pb$subregion, levels = c("LA", "BA", "aBA"))

# =========
# Dream DE
# =========
form = ~ 0 + subregion + (1|subject)

L <- makeContrastsDream(form, colData(pb),
                        contrasts = c(
                        aBA_vs_BA = "subregionaBA - subregionBA",
                        BA_vs_LA = "subregionBA - subregionLA",
                        aBA_vs_LA = "subregionaBA - subregionLA")
                    )

fit <- dream(cobj, form, colData(pb), L)
fit <- eBayes(fit)

fit
#             (Intercept) subregionBA subregionaBA
# CALCR_PENK    -5.287273   -4.968893    -3.145135
# CCK_CNR1      -2.732261    5.357846     3.940437
# LAMP5_EGFR    20.725108   -5.568051    -6.749312
# LAMP5_NOS1     2.889530   -6.366829    -6.354154
# PVALB_MYO5B   -7.214574    2.254768    -4.945953
# 10 more rows ...


pca <- prcomp(t(standardize(cobj)))

# merge with metadata
df_pca <- merge(pca$x, colData(pb), by = "row.names")

# Plot PCA
#   color by Subject
#   shape by Stimulated vs unstimulated
pdf(here(plot_dir, "PCA_subregion.pdf"), width=8, height=6)
ggplot(df_pca, aes(PC1, PC2, color = as.character(subregion))) +
  geom_point(size = 3) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  scale_color_discrete(name = "Subject") +
  xlab("PC1") +
  ylab("PC2")
dev.off()



# ==== build tree =====
hc <- buildClusterTreeFromPB(pb)
res <- treeTest(fit, cobj, hc, coef = "subregionBA")

form <- ~ (1| subregion) + (1|subject) + (1|sex) 
vp <- fitExtractVarPartModel(cobj, form, colData(pb))

library(RColorBrewer)
cols = c(brewer.pal(ncol(vp)-1, "Set1"), "grey85")

p1 <- plotTreeTestBeta(res)
row_order <- ggtree::get_taxa_name(p1)


pdf(here(plot_dir, "forest_BLAsubregions.pdf"), width=3, height=4)
p_forest <- crumblr::plotForest(res, hide = FALSE)  + theme(legend.position = "none") +
  ggtitle("LA <---> BA") +
  theme(axis.text.y = element_blank()) +
  # increase font size for x ticks
  theme(axis.text.x = element_text(size = 14)) 
print(p_forest)
dev.off()


png(here(plot_dir, "forest_BLAsubregion.png"), width=12, height=6)
plotTreeTestBeta(res) +
    theme(legend.position = "bottom", legend.box = "vertical") |
    crumblr::plotForest(res, hide = FALSE) |
    plotPercentBars(vp[row_order,], col=cols) 
dev.off()


# ===========================================================================
# New plotTreeTestBeta function to NOT show tree labels
# ===========================================================================
library(dplyr)
library(ggtree)
plotTreeTestBeta_new <- function(tree, low = "blue", mid = "white", high = "red", xmax.scale = 1.5) {
  # PASS R check
  isTip <- label <- node <- FDR <- NULL

  # comparison only works for fixed effects models
  if (!all(tree@data$method %in% c("FE", "FE.empirical"))) {
    stop("tree1 must be evaluated with a fixed effect model")
  }

  beta_max <- tree %>%
    as_tibble() %>%
    pull(beta) %>%
    abs() %>%
    max()

  fig <- ggtree::ggtree(tree, branch.length = "none") +
    #geom_tiplab(color = "black", size = 4, hjust = 0, offset = .4) + 
    geom_point2(aes(label = node, color = beta, size = pmin(4, -log10(FDR)))) +
    scale_color_gradient2(
      name = bquote(beta), low = low, mid = mid, high = high,
      midpoint = 0, limits = c(-beta_max, beta_max)
    ) +
    scale_size_area(name = bquote(-log[10] ~ FDR), limits = c(0, 4)) +
    geom_text2(aes(label = "*", subset = FDR < 0.05), color = "black", size = 6, vjust = 0.75, hjust = .5) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

  xmax <- layer_scales(fig)$x$range$range[2]

  fig + xlim(0, xmax * xmax.scale)
}


# ===========================================================================
# End new plotTreeTestBeta function
# ===========================================================================




# ================================================
# Heatmap of all contrasts
# ================================================

coef_array = c("BA_vs_LA", "aBA_vs_LA", "aBA_vs_BA")

tab = lapply(coef_array, function(coef){
  
  tab = topTable(fit.contrast, coef=coef, number=Inf)
  tab$coef = coef
  tab$celltype = rownames(tab)
  tab$se = with(tab, logFC/t)

  tab
})


tab = do.call(rbind, tab)
rownames(tab) = c()
tab$adj.P.Val = p.adjust(tab$P.Value, "fdr")

# get order of cell types
res = treeTest( fit.contrast, cobj, hc, coef=coef_array[1])
fig1 = plotTreeTest(res) + theme(legend.position="none") + ggtitle(coef)
lvls = rev(get_taxa_name(fig1))
tab$celltype = factor(tab$celltype, lvls)


# Clean up coefficient labels
tab$coef <- gsub("^subregion", "", tab$coef)
coef_array <- gsub("^subregion", "", coef_array)

# After cleaning coefficient names

# Reassign coef as a factor with correct labels
tab$coef <- factor(gsub("^Source", "", tab$coef), gsub("^Source", "", coef_array))

# Set heatmap limits and aspect ratio
lim <- max(abs(tab$logFC))
ratio <- length(unique(tab$celltype)) / length(unique(tab$coef))

# Plot heatmap
pdf(here(plot_dir, "heatmap_subregion.pdf"), width=6.5, height=5)
p1_heatmap <- ggplot(tab, aes(coef, celltype, fill=logFC, label=ifelse(adj.P.Val < 0.05, "*", ''))) +
  geom_tile() +
  geom_text(vjust=0.75, hjust=0.5, color="black", size=6) +
  theme_classic() + 
  theme(aspect.ratio=ratio, axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient2(low="blue", mid="white", high="red", limits=c(-lim, lim)) +
  ylab(NULL)  +
  theme(axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 11)) +
  theme(legend.position = "left") 

print(p1_heatmap)
dev.off()



# ==============================
# Combined beta and forest plot
# ==============================
pdf(here(plot_dir, "Beta_forest_BLAsubregions.pdf"), width=6, height=5)
    p_beta <- plotTreeTestBeta_new(res, xmax.scale = 1.1) +
        theme(legend.direction = "vertical", legend.box = "vertical") 

    p_forest <- crumblr::plotForest(res, hide = FALSE)  + theme(legend.position = "none") +
    ggtitle("All <---> CN") 

    heatmap <- p1_heatmap + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") 

    print(p_beta + p_forest + heatmap + plot_layout(ncol=3, guides = "collect") + 
    plot_annotation(title = "Cell type composition in CeA vs BLA"))
dev.off()




# ===================================
# Plotting raw proportion 
# ===================================


library(rstatix)
library(dplyr)
library(ggpubr)

# make reguon colors ,"Central Nucleus" = "#FFD700" and BLA = grey
subregion_colors <- c( "Accessory Basal" ="#8A2BE2", "Basal" = "#FF4500", "Lateral" = "#00BFFF")


# Extract cell types, subregions, and sample information
cell_types <- sce$fine_celltype
subregions <- sce$subregion  
samples <- sce$Sample

# Create a data frame with counts per sample
cell_data <- data.frame(
  CellType = cell_types,
  Subregion = subregions,
  Sample = samples
)

# Calculate the proportion of each cell type per sample
cell_proportions <- cell_data %>%
  group_by(Sample, Subregion, CellType) %>%
  summarize(Count = n(), .groups = "drop") %>%  # Add .groups = "drop"
  group_by(Sample, Subregion) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Set the levels of Subregion factor in the desired order
cell_proportions$Subregion <- factor(cell_proportions$Subregion, 
                                     levels = c("Lateral", "Basal", "Accessory Basal", "Central Nucleus"))

# Perform post-hoc Tukey test for each cell type and store results
tukey_results <- cell_proportions %>%
  group_by(CellType) %>%
  tukey_hsd(Proportion ~ Subregion) %>%
  add_significance() %>%
  filter(p.adj.signif != "ns") %>%
  add_xy_position(x="subregion")  # Filter to keep only significant results

# Modify the y.position for plotting (for manual adjustment of p-value locations)
tukey_results <- tukey_results 

# Now plot using ggplot2 and manually add the significant p-values
p3 <- ggplot(cell_proportions, aes(x = Subregion, y = Proportion)) +
  geom_boxplot(aes(fill=Subregion), width = 0.3, outlier.shape = NA) +
  geom_beeswarm(cex = 3, size = 3) +
  #theme_minimal() +
  labs(title = "Proportion of Cell Types per Sample by subregion in Macaques", 
       y = "Proportion", 
       x = "subregion") +
  scale_fill_manual(name = "subregion", values = subregion_colors) +
  facet_wrap(~ CellType, scales = "free") +
      theme_bw()  +
  theme(axis.text.x = #vertical 
        element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 18, hjust = 0.5),
        strip.text = element_text(size = 14),
        # Adjust legend size
        legend.text = element_text(size = 14),  # Increase legend text size
        legend.title = element_text(size = 16),  # Increase legend title size
        legend.key.size = unit(1.5, "lines")) +  # Increase legend key size) +
  #stat_pvalue_manual(tukey_results, label.size=6)


# Save the plot with the significant p-values only
png(here(plot_dir, "Faceted_Boxplot_Excit_Celltype_subregion_Proportion_SignificantOnly.png"), width = 13, height = 13, units = "in", res = 300)
p3
dev.off()



# ==== CGE cell types ====

# Choose both cell types now
celltypes <- c("VIP_SEMA5A", "VIP_NRXN1", "LAMP5_EGFR", "LAMP5_NOS1", "CCK_CNR1", "CALCR_PENK")

# subset cell_proportions
cell_proportions_CeA <- cell_proportions %>%
  filter(CellType %in% celltypes)

# Plot: faceted by subregion (rows) and celltype (columns)
pdf(here(plot_dir, "crumblr_CGE_celltype_proportion_faceted.pdf"), width = 6, height = 3.5)
p3_boxplot_CeA <- ggplot(cell_proportions_CeA, aes(Subregion, Proportion, fill=Subregion)) +
  geom_boxplot(width = .5, outlier.shape = NA) +
  geom_beeswarm(cex = 3, size = 1.5) +
  facet_wrap(~CellType, scales = "free", nrow=2) +
  theme_bw() +
  scale_fill_manual(name = "Subregion", values = subregion_colors) +
  ylab("Proportion") +
  xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("CGE Cell Types in BLA")


print(p3_boxplot_CeA)
dev.off()


# ==== MGE cell types ====

# Choose both cell types now
celltypes <- c("SST_NOS1", "SST_TMEM132C", "SST_NXPH2", "PVALB_MYO5B", "PVALB_UNC5B", "PVALB_ST18")

# subset cell_proportions
cell_proportions_BA <- cell_proportions %>%
  filter(CellType %in% celltypes)

# Plot: faceted by subregion (rows) and celltype (columns)
pdf(here(plot_dir, "crumblr_MGE_celltype_proportion_faceted.pdf"), width = 6, height = 3.5)
p3_boxplot_BA <- ggplot(cell_proportions_BA, aes(Subregion, Proportion, fill=Subregion)) +
  geom_boxplot(width = .5, outlier.shape = NA) +
  geom_beeswarm(cex = 3, size = 1.5) +
  facet_wrap(~CellType, scales = "free", nrow=2) +
  theme_bw() +
  scale_fill_manual(name = "Subregion", values = subregion_colors) +
  ylab("Proportion") +
  xlab(NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
    ggtitle("MGE Cell Types in BLA")
print(p3_boxplot_BA)
dev.off()