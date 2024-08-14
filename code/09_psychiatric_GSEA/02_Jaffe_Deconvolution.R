library("spatialLIBD")
library("patchwork")
library("scran")
library("here")
library("scater")
library("dplyr")
library("tidyr")

# Save directories
plot_dir = here("plots", "09_psychiatric_GSEA", "deconvolution")
processed_dir <- here("processed-data", "09_psychiatric_GSEA")

# ===== Load human marker genes =====
all_markers <- read.csv(here(processed_dir, "top_100_human_markers_fine_celltype.csv"))


# ====== Load human SCE ======
sce <- readRDS(here("processed-data", "09_psychiatric_GSEA", "sce.human_all_genes.rds"))

# filter cell types with < 70 cells
celltype_counts <- table(colData(sce)$fine_celltype)
celltypes_to_drop <- names(celltype_counts[celltype_counts < 70])
sce <- sce[, !(colData(sce)$fine_celltype %in% celltypes_to_drop)]
sce$fine_celltype <- factor(sce$fine_celltype)


# ======= Loading Jaffe bulk data ========

# load Jaffe bulk data
load(here(processed_dir, "rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata"))

# subset to basolateral amygdala
rownames(rse_gene) <- rowData(rse_gene)$Symbol
rse_amy <- rse_gene[,rse_gene$Region == "BasoAmyg"]
dim(rse_amy)

# subset to MDD, PTSD, NTC
rse_amy <- rse_amy[, rse_amy$PrimaryDx %in% c("PTSD","MDD", "Control")]
rse_amy <- rse_amy[, rse_amy$PrimaryDx != "NA"]
unique(rse_amy$PrimaryDx)

# reset factors
rse_amy$PrimaryDx <- factor(rse_amy$PrimaryDx)
rse_amy
# class: RangedSummarizedExperiment 
# dim: 58037 294 
# metadata(0):
# assays(1): counts
# rownames(58037): DDX11L1 WASH7P ... MT-TT MT-TP
# rowData names(10): Length gencodeID ... NumTx gencodeTx
# colnames(294): R16254 R15943 ... R17169 R17173
# colData names(70): RNum BrNum ... gene_Unassigned_Duplicate rRNA_rate


# ========= Deconvolution ==========

unique(colData(sce)$fine_celltype)
#  [1] LAMP5_NTNG1     ESR1_ADRA1A     VIP_ADRA1B      MEIS2_COL25A1  
#  [5] THSD7B_CALB2    MEIS1_PARD3B    TSHZ1.1         ZBTB20_SLC4A4  
#  [9] SST_NOS1        PEX5L_MYRIP     PVALB_MYO5B     Astrocyte      
# [13] SST_PRKCQ       LAMP5_KIT       Oligodendrocyte Microglia      
# [17] CCK_CNR1        ZFHX3_SCN5A     PVALB_UNC5B     Endothelial    
# [21] ST18_IL1RAPL2   LHX8_ANGPT1     ST18_ABCA8      VIP_PLPP4      
# [25] TSHZ1.2         ADARB2_TRPS1    GULP1_TRHDE     RXFP1_KIAA1217 
# [29] Ependymal       SATB2_MPPED1   



# ========= Bulk deconvolution of broad celltypes =========

# get top 25 markers genes for each cell type
marker_genes <- all_markers |>
    filter(std.logFC.rank<= 25, gene %in% rownames(rse_amy)) |>
    pull(gene)

marker_genes <- unique(marker_genes)


# prep data
exp_set_bulk <- Biobase::ExpressionSet(
    assayData = assays(rse_amy[marker_genes, ])$counts,
    phenoData = AnnotatedDataFrame(
        as.data.frame(colData(rse_amy))[c("SAMPLE_ID")]
    )
)

exp_set_sce <- Biobase::ExpressionSet(
    assayData = as.matrix(assays(sce[marker_genes, ])$counts),
    phenoData = AnnotatedDataFrame(
        as.data.frame(colData(sce))[, c("fine_celltype", "Sample")]
    )
)


## check for nuclei with 0 marker expression
zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ", sum(!zero_cell_filter), " cells")
# Exclude 19 cells


exp_set_sce <- exp_set_sce[, zero_cell_filter]

# Reference-based deconvolution
est_prop <- BisqueRNA::ReferenceBasedDecomposition(
    bulk.eset = exp_set_bulk,
    sc.eset = exp_set_sce,
    cell.types = "fine_celltype",
    subject.names = "Sample",
    use.overlap = FALSE
)

# exploring outputs
est_prop$bulk.props <- t(est_prop$bulk.props)

## add metadata to proportion estimates
pd <- colData(rse_amy) |>
    as.data.frame() |>
    select(Sample = RNum, Sex, PrimaryDx)

## make proption estimates long so they are ggplot friendly
prop_long <- est_prop$bulk.props |>
    as.data.frame() |>
    tibble::rownames_to_column("Sample") |>
    tidyr::pivot_longer(!Sample, names_to = "cell_type", values_to = "prop") |>
    left_join(pd)
prop_long
# # A tibble: 6 × 5
#   Sample cell_type       prop Sex   PrimaryDx
#   <chr>  <chr>          <dbl> <chr> <fct>    
# 1 R16254 ADARB2_TRPS1 0.0568  M     PTSD     
# 2 R16254 Astrocyte    0.0604  M     PTSD     
# 3 R16254 CCK_CNR1     0.00707 M     PTSD     
# 4 R16254 Endothelial  0.0145  M     PTSD     
# 5 R16254 Ependymal    0.0229  M     PTSD     
# 6 R16254 ESR1_ADRA1A  0.0538  M     PTSD    

# dropp NA from PrimaryDx
prop_long <- prop_long[!is.na(prop_long$PrimaryDx),]


png(here(plot_dir, "bulk_deconvolution_boxplot.png"), width = 10, height = 10, units = "in", res = 300)
ggplot(prop_long, aes(x = cell_type, y = prop, fill = PrimaryDx)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "Cell Type", y = "Proportion", fill = "Diagnosis") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()





## create composition bar plots
png(here(plot_dir, "bulk_deconvolution_composition_bar.png"), width = 10, height = 10, units = "in", res = 300)
plot_composition_bar(prop_long = prop_long, sample_col = "Sample", x_col = "PrimaryDx")
dev.off()
# # A tibble: 8,820 × 5
#    Sample cell_type       prop Sex   PrimaryDx
#    <chr>  <chr>          <dbl> <chr> <fct>    
#  1 R16254 ADARB2_TRPS1 0.0568  M     PTSD     
#  2 R16254 Astrocyte    0.0604  M     PTSD     
#  3 R16254 CCK_CNR1     0.00707 M     PTSD     
#  4 R16254 Endothelial  0.0145  M     PTSD     
#  5 R16254 Ependymal    0.0229  M     PTSD     
#  6 R16254 ESR1_ADRA1A  0.0538  M     PTSD     
#  7 R16254 GULP1_TRHDE  0.185   M     PTSD     
#  8 R16254 LAMP5_KIT    0.0325  M     PTSD     
#  9 R16254 LAMP5_NTNG1  0.0294  M     PTSD     
# 10 R16254 LHX8_ANGPT1  0.0247  M     PTSD  
