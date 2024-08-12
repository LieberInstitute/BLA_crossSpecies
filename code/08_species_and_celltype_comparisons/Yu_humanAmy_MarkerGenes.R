library(Seurat)
library(SingleCellExperiment)
library(scran)
library(here)
library(DeconvoBuddies)
library(dplyr)

# directories
processed_dir = here("processed-data")
plot_dir = here("plots", "08_species_comparisons")

yu.seurat <- readRDS(here(processed_dir, "Yu_humanAmy_seurat.rds"))
yu.seurat
# An object of class Seurat 
# 33538 features across 91699 samples within 1 assay 
# Active assay: RNA (33538 features, 0 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

# convert to SCE
yu.sce <- as.SingleCellExperiment(yu.seurat)
yu.sce
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


colnames(colData(yu.sce))
# [1] "barcode"            "library_id"         "sample"             "n_genes_by_counts"  "total_counts"      
# [6] "total_counts_mt"    "pct_counts_mt"      "doublet_scores"     "predicted_doublets" "celltype"          
# [11] "clusters"           "space_anno"         "orig_anno"          "ident"  

unique(yu.sce$ident)
# [1] Human_Endo NOSTRIN   Human_LAMP5 ABO      Human_Astro_1 FGFR3  Human_PVALB ADAMTS5  Human_LAMP5 COL25A1 
# [6] Human_SOX11 EBF2     Human_Oligo_5 OPALIN Human_VGLL3 CNGB1    Human_LAMP5 COL14A1  Human_Oligo_1 OPALIN
# [11] Human_Astro_4 FGFR3  Human_Astro_3 FGFR3  Human_CALCR LHX8     Human_HGF ESR1       Human_HGF C11orf87  
# [16] Human_HGF NPSR1      Human_DRD2 PAX6      Human_VGLL3 MEPE     Human_PRKCD          Human_TSHZ1 CALCRL  
# [21] Human_LAMP5 BDNF     Human_VIP ABI3BP     Human_HTR3A DRD2     Human_VIP NDNF       Human_LAMP5 NDNF    
# [26] Human_Astro_2 FGFR3  Human_OPC_3 PDGFRA   Human_DRD2 ISL1      Human_SST HGF        Human_SATB2 CALCRL  
# [31] Human_TFAP2C         Human_SATB2 IL15     Human_Oligo_2 OPALIN Human_OPC_2 PDGFRA   Human_TSHZ1 SEMA3C  
# [36] Human_Micro CTSS     Human_SST EPYC       Human_RXFP2 RSPO2    Human_Oligo_3 OPALIN Human_OPC_4 PDGFRA  
# [41] Human_SATB2 ST8SIA2  Human_STRIP2         Human_OPC_1 PDGFRA   Human_Oligo_4 OPALIN Human_Oligo_6 OPALIN
# 45 Levels: Human_Oligo_5 OPALIN Human_Oligo_6 OPALIN Human_Astro_2 FGFR3 Human_LAMP5 COL25A1 ... Human_OPC_1 PDGFRA

unique(yu.sce$celltype)
# [1] Endothelial     ExN             Astrocyte       InN             Oligodendrocyte OPC             Microglia      
# Levels: Astrocyte Endothelial ExN InN Microglia OPC Oligodendrocyte

# ====================
#  Marker Genes
# ====================


# === Broad celltypes ===

# 1 vs All
oneVall <- findMarkers_1vAll(
    yu.sce,
    assay_name = "logcounts",
    cellType_col = "celltype",
    add_symbol = FALSE,
    mod = "~sample",
    verbose = TRUE
)

top_markers <- oneVall %>%
    group_by(cellType.target) %>%
    top_n(n = 100, wt = anno_logFC)

write.csv(top_markers, here(processed_dir, "top_100_markers_broad_celltype.csv"), row.names = FALSE)

# mean ratio markers
ratio_markers <- get_mean_ratio2(
    yu.sce,
    cellType_col = "celltype",
    assay_name = "logcounts",
    add_symbol = TRUE
)

top_ratio_markers <- ratio_markers %>%
    group_by(cellType.target) %>%
    top_n(n = 100, wt = ratio)

write.csv(top_ratio_markers, here(processed_dir, "top_100_ratio_markers_broad_celltype.csv"), row.names = FALSE)


# === Fine celltypes ===

# 1 vs All
oneVall <- findMarkers_1vAll(
    yu.sce,
    assay_name = "logcounts",
    cellType_col = "ident",
    add_symbol = FALSE,
    mod = "~sample",
    verbose = TRUE
)

top_markers <- oneVall %>%
    group_by(cellType.target) %>%
    top_n(n = 100, wt = anno_logFC)

write.csv(top_markers, here(processed_dir, "top_100_markers_fine_celltype.csv"), row.names = FALSE)

# mean ratio markers
ratio_markers <- get_mean_ratio2(
    yu.sce,
    cellType_col = "ident",
    assay_name = "logcounts",
    add_symbol = TRUE
)

top_ratio_markers <- ratio_markers %>%
    group_by(cellType.target) %>%
    top_n(n = 100, wt = ratio)

write.csv(top_ratio_markers, here(processed_dir, "top_100_ratio_markers_fine_celltype.csv"), row.names = FALSE)
