library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")
library("scuttle")
library("SpatialExperiment")
library("here")



#Load the object
load("sce_FINAL_all_celltypes.rda", verbose = TRUE)
#sce <- rda
# [1] "orig.ident"             "nCount_originalexp"     "nFeature_originalexp"   "Sample_num"            
# [5] "Sample"                 "Species"                "Subject"                "Sex"                   
# [9] "Region"                 "Subregion"              "DV_axis"                "PI.NeuN"               
# [13] "batch"                  "Barcode"                "sum"                    "detected"              
# [17] "subsets_Mito_sum"       "subsets_Mito_detected"  "subsets_Mito_percent"   "total"                 
# [21] "high_mito"              "low_lib"                "low_genes"              "discard_auto"          
# [25] "doubletScore"           "sizeFactor"             "species"                "unintegrated_clusters" 
# [29] "seurat_clusters"        "integrated_snn_res.0.5" "ident"                  "key"                   
# [33] "broad_celltype"         "fine_celltype"   



colors <- readRDS("celltype_colors.rds")
colors
# LAMP5_NTNG1      VIP_ADRA1B    THSD7B_CALB2         TSHZ1.1        SST_NOS1     PVALB_MYO5B 
# "#1F78C8"       "#ff0000"       "#33a02c"       "#6A33C2"       "#ff7f00"       "#565656" 
# PRKCD_DRD2       SST_PRKCQ       LAMP5_KIT        CCK_CNR1     ZFHX3_SCN5A     PVALB_UNC5B 
# "#FFD700"       "#a6cee3"       "#FB6496"       "#b2df8a"       "#CAB2D6"       "#FDBF6F" 
# CARTPT_CDH23   ST18_IL1RAPL2     LHX8_ANGPT1       VIP_PLPP4      PRKCD_DRD1         TSHZ1.2 
# "#999999"       "#EEE685"       "#C8308C"       "#FF83FA"       "#C814FA"       "#0000FF" 
# ESR1_ADRA1A   MEIS2_COL25A1    MEIS1_PARD3B   ZBTB20_SLC4A4     PEX5L_MYRIP      ST18_ABCA8 
# "#E41A1C"       "#66628D"       "#419486"       "#5A9D5A"       "#91569A"       "#D96D3B" 
# ADARB2_TRPS1     GULP1_TRHDE  RXFP1_KIAA1217      GRIK3_TNS3    SATB2_MPPED1 SLC17A8_ST8SIA2 
# "#FFAD12"       "#F6EF32"       "#B6742A"       "#D26D7A"       "#DD87B4"       "#999999" 
# Astrocyte Oligodendrocyte       Microglia     Endothelial       Ependymal             OPC 
# "#4D4D4D"       "#7F7F7F"       "#A0A0A0"       "#BBBBBB"       "#D1D1D1"       "#E6E6E6" 

#Change the rownames frome ensembl id to gene_name
#rownames(sce) <-  uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)

#Source
source("initial.R", print.eval = TRUE)

#increase max number of colors
rda <- registerAppOptions(rda, color.maxlevels = length(unique(rda$fine_celltype)))

#Deploy app
iSEE(
  rda,
  appTitle = "BLA CrossSpecies",
  initial = initial,
  colormap = ExperimentColorMap(
      colData = list(fine_celltype = function(x) {
          colors                       # Return color mapping
      })  # Pass the function for fine_celltype colors
  )
)

