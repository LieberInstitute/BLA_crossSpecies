library(dplyr)
library(here)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(tidyr)
library(scater)
library(BisqueRNA)
library(DeconvoBuddies)
library(tidySingleCellExperiment)
library(tidySummarizedExperiment)
library(tidybulk)


# directories
processed_dir = here("processed-data", "10_bulk_deconvolution")
plot_dir = here("plots", "10_bulk_deconvolution")

# load data
sce <- readRDS(here("processed-data", "Yu_humanAmy_seurat.rds"))
sce <- Seurat::as.SingleCellExperiment(sce)
sce
# class: SingleCellExperiment 
# dim: 33538 91699 
# metadata(0):
#     assays(2): counts logcounts
# rownames(33538): MIR1302-2HG FAM138A ... AC213203.1 FAM231C
# rowData names(0):
#     colnames(91699): AAACCCAAGTCGCCAC-1 AAACCCAGTCGAGATG-1 ... TTTGTTGTCACTTCTA-14 TTTGTTGTCCAAATGC-14
# colData names(14): barcode library_id ... orig_anno ident
# reducedDimNames(2): PCA UMAP
# mainExpName: RNA
# altExpNames(0):

load(here("processed-data", "rse_gene_PTSD_VA_LIBD_qcAndAnnotated_n1285.Rdata"))
rse_gene
# class: RangedSummarizedExperiment 
# dim: 58037 1285 
# metadata(0):
#     assays(1): counts
# rownames(58037): ENSG00000223972.5 ENSG00000227232.5 ... ENSG00000210195.2 ENSG00000210196.2
# rowData names(10): Length gencodeID ... NumTx gencodeTx
# colnames(1285): R16254 R16050 ... R17174 R17175
# colData names(70): RNum BrNum ... gene_Unassigned_Duplicate rRNA_rate


#  =========== subsetting to amygdala sample =========
rownames(rse_gene) <- rowData(rse_gene)$Symbol
colnames(colData(rse_gene))
# [1] "RNum"                           "BrNum"                          "RIN"                            "Region"                        
# [5] "AgeDeath"                       "Sex"                            "Race"                           "PrimaryDx"                     
# [9] "Group"                          "plate"                          "Position"                       "SAMPLE_ID"                     
# [13] "FQCbasicStats"                  "perBaseQual"                    "perTileQual"                    "perSeqQual"                    
# [17] "perBaseContent"                 "GCcontent"                      "Ncontent"                       "SeqLengthDist"                 
# [21] "SeqDuplication"                 "OverrepSeqs"                    "AdapterContent"                 "KmerContent"                   
# [25] "SeqLength_R1"                   "percentGC_R1"                   "phred20-21_R1"                  "phred48-49_R1"                 
# [29] "phred76-77_R1"                  "phred100-101_R1"                "phredGT30_R1"                   "phredGT35_R1"                  
# [33] "Adapter50-51_R1"                "Adapter70-71_R1"                "Adapter88-89_R1"                "SeqLength_R2"                  
# [37] "percentGC_R2"                   "phred20-21_R2"                  "phred48-49_R2"                  "phred76-77_R2"                 
# [41] "phred100-101_R2"                "phredGT30_R2"                   "phredGT35_R2"                   "Adapter50-51_R2"               
# [45] "Adapter70-71_R2"                "Adapter88-89_R2"                "ERCCsumLogErr"                  "bamFile"                       
# [49] "trimmed"                        "numReads"                       "numMapped"                      "numUnmapped"                   
# [53] "overallMapRate"                 "concordMapRate"                 "totalMapped"                    "mitoMapped"                    
# [57] "mitoRate"                       "totalAssignedGene"              "gene_Assigned"                  "gene_Unassigned_Ambiguity"     
# [61] "gene_Unassigned_MultiMapping"   "gene_Unassigned_NoFeatures"     "gene_Unassigned_Unmapped"       "gene_Unassigned_MappingQuality"
# [65] "gene_Unassigned_FragmentLength" "gene_Unassigned_Chimera"        "gene_Unassigned_Secondary"      "gene_Unassigned_Nonjunction"   
# [69] "gene_Unassigned_Duplicate"      "rRNA_rate" 

unique(rse_gene$Region)
# [1] "BasoAmyg"   "MedialAmyg" "dACC"       "DLPFC"  

rse_amy <- rse_gene[,rse_gene$Region == "BasoAmyg"]
dim(rse_amy)
# [1] 58037   322

rse_amy <- rse_amy[, rse_amy$PrimaryDx %in% c("PTSD","MDD", "Control")]
rse_amy <- rse_amy[, rse_amy$PrimaryDx != "NA"]
unique(rse_amy$PrimaryDx)
# [1] "PTSD"    "Control" "MDD"  

# reset factors
rse_amy$PrimaryDx <- factor(rse_amy$PrimaryDx)

# ============ Marker genes of Yu et al ============

unique(sce$celltype)
# [1] Endothelial     ExN             Astrocyte       InN             Oligodendrocyte OPC             Microglia   

# load in marker genes from .csv files
ratio_markers <- read.csv(here("processed-data", 
                              "08_species_comparisons", 
                              "Yu_humanAmy_MarkerGenes", 
                              "top_100_ratio_markers_broad_celltype.csv")
                         )

all_markers <- read.csv(here("processed-data", 
                             "08_species_comparisons", 
                             "Yu_humanAmy_MarkerGenes", 
                             "top_100_markers_broad_celltype.csv")
                        )


# === Broad markers ===
# plot broad markers using ratio
DeconvoBuddies::plot_marker_express(sce,
                                    stats = all_markers,
                                    cell_type = "InN",
                                    cellType_col = "celltype",
                                    n_genes = 10,
                                    rank_col = "rank_marker",
                                    anno_col = "anno_logFC",
) + 
    ggtitle("InN Top 10 Markers: 1vAll")

# plot broad markers using 1vAll
DeconvoBuddies::plot_marker_express(sce,
                                    stats = ratio_markers,
                                    cell_type = "InN",
                                    cellType_col = "celltype",
                                    n_genes = 10,
                                    rank_col = "rank_ratio",
                                    anno_col = "anno_ratio",
) + 
    ggtitle("InN Top 10 Markers: Ratio")




# ========= Bulk deconvolution of broad celltypes =========

# get top 25 markers genes for each cell type
marker_genes <- ratio_markers |>
    filter(rank_ratio <= 5, gene %in% rownames(rse_amy)) |>
    pull(gene)

marker_genes <- unique(marker_genes)

# 
# plotGroupedHeatmap(sce, features=marker_genes,group="celltype",block="sample", center=TRUE, scale=TRUE)
# 
# sce.pb <- aggregateAcrossCells(sce, ids = colData(sce)[,c("celltype", "sample")], use.assay.type = "logcounts")
# sce.pb$pb_samples <- paste0(sce.pb$sample, "_", sce.pb$celltype)
# 
# plotHeatmap(sce.pb, features=marker_genes, center=TRUE, scale=TRUE)

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
        as.data.frame(colData(sce))[, c("celltype", "sample")]
    )
)

## check for nuclei with 0 marker expression
zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ", sum(!zero_cell_filter), " cells")
#> Exclude 0 cells

exp_set_sce <- exp_set_sce[, zero_cell_filter]

# Reference-based deconvolution
est_prop <- ReferenceBasedDecomposition(
    bulk.eset = exp_set_bulk,
    sc.eset = exp_set_sce,
    cell.types = "celltype",
    subject.names = "sample",
    use.overlap = FALSE
)

# exploring outputs
est_prop$bulk.props <- t(est_prop$bulk.props)
head(est_prop$bulk.props)
#         Astrocyte Endothelial       ExN       InN  Microglia        OPC Oligodendrocyte
# R16254 0.07573618  0.00000000 0.3722180 0.2835025 0.00000000 0.02024928      0.24829399
# R15943 0.11509558  0.02197533 0.3750479 0.3228760 0.02417396 0.04245070      0.09838054
# R15944 0.06053986  0.00000000 0.3724725 0.2960113 0.00000000 0.00484034      0.26613603
# R15945 0.04400988  0.00000000 0.2381548 0.1671616 0.00000000 0.01721521      0.53345853
# R15946 0.07841070  0.01283352 0.3162870 0.2892983 0.01030455 0.03967416      0.25319175
# R15947 0.07463561  0.01522881 0.4040842 0.2552598 0.01733160 0.02721363      0.20624635

## add Phenotype data to proportion estimates
pd <- colData(rse_amy) |>
    as.data.frame() |>
    select(Sample = RNum, Sex, PrimaryDx)

## make proption estimates long so they are ggplot friendly
prop_long <- est_prop$bulk.props |>
    as.data.frame() |>
    tibble::rownames_to_column("Sample") |>
    tidyr::pivot_longer(!Sample, names_to = "cell_type", values_to = "prop") |>
    left_join(pd)

## create composition bar plots
plot_composition_bar(prop_long = prop_long, sample_col = "Sample", x_col = "PrimaryDx")



# =========== Fine cell-types =============


unique(sce$celltype)
# [1] Endothelial     ExN             Astrocyte       InN             Oligodendrocyte OPC             Microglia   

# load in marker genes from .csv files
ratio_markers <- read.csv(here("processed-data", 
                               "08_species_comparisons", 
                               "Yu_humanAmy_MarkerGenes", 
                               "top_100_ratio_markers_fine_celltype.csv")
)

all_markers <- read.csv(here("processed-data", 
                             "08_species_comparisons", 
                             "Yu_humanAmy_MarkerGenes", 
                             "top_100_markers_fine_celltype.csv")
)

#drop any genes startign with AC0
ratio_markers <- ratio_markers[!grepl("^AC0", ratio_markers$gene),]
all_markers <- all_markers[!grepl("^AC0", all_markers$gene),]

# === Broad markers ===
# plot broad markers using ratio
DeconvoBuddies::plot_marker_express(sce,
                                    stats = all_markers,
                                    cell_type = "Human_SST EPYC",
                                    cellType_col = "ident",
                                    n_genes = 5,
                                    rank_col = "rank_marker",
                                    anno_col = "anno_logFC",
) + 
    ggtitle("SST_HGF Top 5 Markers: 1vAll")

# plot broad markers using 1vAll
DeconvoBuddies::plot_marker_express(sce,
                                    stats = ratio_markers,
                                    cell_type = "Human_SST EPYC",
                                    cellType_col = "ident",
                                    n_genes = 5,
                                    rank_col = "rank_ratio",
                                    anno_col = "anno_ratio",
) + 
    ggtitle("SST_HGF Top 5 Markers: Ratio")


# ========= Bulk deconvolution of broad celltypes =========

# get top 25 markers genes for each cell type
marker_genes <- all_markers |>
    filter(rank_marker<= 5, gene %in% rownames(rse_amy)) |>
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
        as.data.frame(colData(sce))[, c("ident", "sample")]
    )
)

## check for nuclei with 0 marker expression
zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ", sum(!zero_cell_filter), " cells")
#> Exclude 0 cells

exp_set_sce <- exp_set_sce[, zero_cell_filter]

# Reference-based deconvolution
est_prop <- ReferenceBasedDecomposition(
    bulk.eset = exp_set_bulk,
    sc.eset = exp_set_sce,
    cell.types = "ident",
    subject.names = "sample",
    use.overlap = FALSE
)


# exploring outputs
est_prop$bulk.props <- t(est_prop$bulk.props)
head(est_prop$bulk.props)
#         Astrocyte Endothelial       ExN       InN  Microglia        OPC Oligodendrocyte
# R16254 0.07573618  0.00000000 0.3722180 0.2835025 0.00000000 0.02024928      0.24829399
# R15943 0.11509558  0.02197533 0.3750479 0.3228760 0.02417396 0.04245070      0.09838054
# R15944 0.06053986  0.00000000 0.3724725 0.2960113 0.00000000 0.00484034      0.26613603
# R15945 0.04400988  0.00000000 0.2381548 0.1671616 0.00000000 0.01721521      0.53345853
# R15946 0.07841070  0.01283352 0.3162870 0.2892983 0.01030455 0.03967416      0.25319175
# R15947 0.07463561  0.01522881 0.4040842 0.2552598 0.01733160 0.02721363      0.20624635

## add Phenotype data to proportion estimates
pd <- colData(rse_amy) |>
    as.data.frame() |>
    select(Sample = RNum, Sex, PrimaryDx)

## make proption estimates long so they are ggplot friendly
prop_long <- est_prop$bulk.props |>
    as.data.frame() |>
    tibble::rownames_to_column("Sample") |>
    tidyr::pivot_longer(!Sample, names_to = "cell_type", values_to = "prop") |>
    left_join(pd)

## create composition bar plots
plot_composition_bar(prop_long = prop_long, sample_col = "Sample", x_col = "PrimaryDx")
prop_long
# Sample cell_type              prop Sex   PrimaryDx
# <chr>  <chr>                 <dbl> <chr> <fct>    
#     1 R16254 Human_Oligo_5 OPALIN 0.0114 M     PTSD     
# 2 R16254 Human_Oligo_6 OPALIN 0.121  M     PTSD     
# 3 R16254 Human_Astro_2 FGFR3  0.0209 M     PTSD     
# 4 R16254 Human_LAMP5 COL25A1  0.0954 M     PTSD     
# 5 R16254 Human_Oligo_4 OPALIN 0.0720 M     PTSD     
# 6 R16254 Human_LAMP5 BDNF     0      M     PTSD     
# 7 R16254 Human_OPC_3 PDGFRA   0.0389 M     PTSD     
# 8 R16254 Human_LAMP5 NDNF     0.0223 M     PTSD     
# 9 R16254 Human_VGLL3 MEPE     0.0213 M     PTSD     
# 10 R16254 Human_HGF NPSR1      0.0922 M     PTSD   

# do fold change instead
prop_diff <- prop_long |>
    group_by(cell_type, PrimaryDx) |>
    summarise(mean_prop = mean(prop), .groups = "drop") |>
    pivot_wider(names_from = PrimaryDx, values_from = mean_prop) |>
    mutate(diff = PTSD - Control) |>
    arrange(desc(diff))

#plot
ggplot(prop_diff, aes(x = reorder(cell_type, diff), y = diff, fill = diff > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Cell Type", y = "log2(PTSD/Control)") +
    ggtitle("PTSD - Control")




assays(rse_amy)$logcounts <- log2(assay(rse_amy) + 1)
genes <- c("SST", "CORT","NPY", "PVALB", "VIP", "CCK")
plot_gene_express(rse_amy, genes, cat="PrimaryDx", assay_name="logcounts")
