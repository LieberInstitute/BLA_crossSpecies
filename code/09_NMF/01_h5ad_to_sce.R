library(biomaRt)
library(zellkonverter)
library(SingleCellExperiment)
library(RcppML)
library(here)

# ==== set up directories ====
processed_dir <- here("processed-data")

# ==== load in retro seq data ====
amy_h5ad <- readH5AD(here(processed_dir,"retro-seq-data","amy_mch_matrix.h5ad"))
orthology <- read.csv(here(processed_dir, "retro-seq-data","human_mouse_orthologs.csv"))

# ======= Converting gene names =======
mart = useMart('ensembl')

# get mouse databse
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# get symbols
symb <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),
              filters = "ensembl_gene_id", 
              values = rownames(amy_h5ad),
              mart = mart)

# get only matches
symbs <- symb$mgi_symbol[match(rownames(amy_h5ad), symb$ensembl_gene_id, nomatch = NA)]

# add gene names to rowData
rowData(amy_h5ad)$gene_name<-symbs

# drop start and end rowData (?)
rowData(amy_h5ad)$start<-NULL
rowData(amy_h5ad)$end<-NULL

# Drop NAs and set rownames as gene_names
amy_h5ad<-amy_h5ad[!is.na(rowData(amy_h5ad)$gene_name),]
rownames(amy_h5ad)<-rowData(amy_h5ad)$gene_name

# check is NA in rownames
sum(is.na(rownames(amy_h5ad)))
#check is NA in orthology column 3
sum(is.na(orthology$Column3))

# ======= Convert to human orthologs =======

# drop any genes in orthology not in AMY
names <- orthology[orthology$Column3 %in% rownames(amy_h5ad),]
names <- names[match(rownames(amy_h5ad), names$Column3),]

# drop excess genes in amy
amy_h5ad <- amy_h5ad[rownames(amy_h5ad) %in% names$Column3,]

# check that the orders are the same
head(rownames(amy_h5ad))
# [1] "Xkr4"   "Rp1"    "Sox17"  "Mrpl15" "Lypla1" "Tcea1" 
head(names$Column3)
# [1] "Xkr4"   "Rp1"    "Sox17"  "Mrpl15" "Lypla1" "Tcea1" 

# add human orthologs to rowData
rowData(amy_h5ad)$human_ortholog<-names$Column1

# save SCE
saveRDS(amy_h5ad, here(processed_dir,"retro-seq-data","amy_sce_with_orthologs.rds"))

sce <- amy_h5ad
sce


# load
sce <- readRDS(here(processed_dir,"retro-seq-data","amy_sce_with_orthologs.rds"))
sce
# class: SingleCellExperiment 
# dim: 19507 2003 
# metadata(1): log
# assays(1): X
# rownames(19507): Xkr4 Rp1 ... mt-Nd6 mt-Cytb
# rowData names(3): chrom gene_name human_ortholog
# colnames(2003): Pool154_Plate1-1-A1-A13 Pool154_Plate1-1-A1-A14 ... Pool180_Plate13-6-I24-P23
# Pool180_Plate13-6-I24-P24
# colData names(0):
#     reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

# read in csv with metadata
metadata <- read.csv(here(processed_dir,"retro-seq-data","cell_48032_RS2_meta_nooutlier.csv"))
head(metadata)
# 1 180118_CEMBA_mm_P56_P63_RS2_180117_P1_180117_P2_A10_AD001 0.006500221 0.7758496 0.03202742      1090870 180117_P1
# 2 180118_CEMBA_mm_P56_P63_RS2_180117_P1_180117_P2_A10_AD002 0.006709048 0.7751456 0.03409142      1458065 180117_P1
# 3 180118_CEMBA_mm_P56_P63_RS2_180117_P1_180117_P2_A10_AD004 0.007130547 0.7758953 0.03589979      1614993 180117_P1
# 4 180118_CEMBA_mm_P56_P63_RS2_180117_P1_180117_P2_A10_AD006 0.006032079 0.7733448 0.02859073      1598778 180117_P1
# 5 180118_CEMBA_mm_P56_P63_RS2_180117_P1_180117_P2_A10_AD007 0.005006534 0.7734715 0.02202152       762007 180117_P2
# 6 180118_CEMBA_mm_P56_P63_RS2_180117_P1_180117_P2_A10_AD008 0.005096293 0.7698228 0.02115723      1872338 180117_P2
# PlateNormReads Experiment Source Target    Sex  L1     L2        L3           L4 L1.annot PassTargetFilter Cocluster
# 1   -0.555116130       Tm3C    MOp     SC   male c14 c14_c6 c14_c6_c0 c14_c6_c0_c0       ET             TRUE    CTX-13
# 2   -0.136540279       Tf3C    MOp     SC female c14 c14_c5 c14_c5_c0 c14_c5_c0_c0       ET             TRUE    CTX-13
# 3    0.010932597       Tm4B    MOp     SC   male c14 c14_c6 c14_c6_c1 c14_c6_c1_c0       ET             TRUE    CTX-13
# 4   -0.003625689       Tf4B    MOp     SC female c14 c14_c0 c14_c0_c0 c14_c0_c0_c1       ET             TRUE    CTX-13
# 5   -1.169393312       Pm3C    MOp    STR   male  c1  c1_c1  c1_c1_c3  c1_c1_c3_c2   IT-Sup             TRUE     CTX-2
# 6    0.127571431       Pf3C    MOp    STR female  c1  c1_c1  c1_c1_c3  c1_c1_c3_c0   IT-Sup             TRUE     CTX-2
# Subclass
# 1   L5 ET CTX Glut
# 2   L5 ET CTX Glut
# 3   L5 ET CTX Glut
# 4   L5 ET CTX Glut
# 5 L4/5 IT CTX Glut
# 6 L4/5 IT CTX Glut

# === Adding metadata ===
amy.metadata <- metadata[metadata$Source == "AMY",]
colData(sce)$Cell.Name <- colnames(sce)
test <- plyr::join(as.data.frame(colData(sce)), amy.metadata, by="Cell.Name", type="inner")
colData(sce) <- DataFrame(test)

# === up object ===

assay(sce, "logcounts") <- assay(sce,"X")
assay(sce, "X") <- NULL
colnames(sce) <- sce$Cell.Name

sce<-sce[,!is.na(sce$Subclass)]

# remove subclass == ""
sce <- sce[,sce$Subclass != ""]


# === Normalization, Dim reduction, and Clustering === 

# Normalization
library(scran)
# set.seed(100)
# clust <- quickCluster(sce) 
# sce.test <- calculateSumFactors(sce, cluster=clust, min.mean=0.1)
# sce <- logNormCounts(sce)

library(harmony)
sce <- RunHarmony(sce, c("Experiment", "Sex"))

# feature selection
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n=4000)

# dim Reduction
sce <- runPCA(sce, ncomponents = 50, subset_row = hvg)
sce <- runUMAP(sce, dimred = "HARMONY", min_dist = 0.5)


plotReducedDim(sce, "UMAP", colour_by = "Target")


# For each target region, get percent of Subclass
data <-data.frame(sce$Subclass,
                          sce$Target)


# get percent of Subclass in each Target
library(dplyr)
data.percent <- data %>% 
    group_by(Target, Subclass) %>% 
    summarise(n=n()) %>% 
    mutate(percent=n/sum(n))

# plot stacked bars for each target
library(ggplot2)
ggplot(data.percent, aes(x=Target, y=percent, fill=Subclass)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="Subclass Percent in Target Regions", x="Target", y="Percent")

# use more distinct colors. I need 27
library(RColorBrewer)
nb.cols <- 27
mycolors <- colorRampPalette(brewer.pal(9, "Set3"))(nb.cols)

ggplot(data.percent, aes(x=Target, y=percent, fill=Subclass)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title="Subclass Percent in Target Regions", x="Target", y="Percent") +
    scale_fill_manual(values = mycolors)
library(bluster)
clust2 <- clusterRows(reducedDim(sce, "HARMONY"), NNGraphParam())
colData(sce)$graph_clust <- factor(clust2)

plotReducedDim(sce, "UMAP", colour_by = "Subclass")

features <- c("Slc17a7", "Gad2", "Cyp26b1","Sst","Rspo2","Meis2")
plotExpression(sce, features, x = "graph_clust", colour_by = "graph_clust")
