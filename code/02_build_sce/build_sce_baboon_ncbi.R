# qrsh -l mem_free=80G,h_vmem=80G

library("SingleCellExperiment")
library("DropletUtils")
library("here")
library("rtracklayer")
library("sessioninfo")
library("dplyr")

# Read in the CSV file as a data frame
tmp <- read.delim(here("raw-data","sampleinfo",
                       "master_sampleinfo_2023-10-09.csv"),
                  header = T,sep=',')

# View the subsetted data
head(tmp)
#        Sample Species Subject Sex   Region Subregion DV_axis  PI.NeuN
#   1 3c-AMYBLA   Human  Br2327     Amygdala       BLA         PI+NeuN+
#   2 4c-AMYBLA   Human  Br8692     Amygdala       BLA         PI+NeuN+
#   3 5c-AMYBLA   Human  Br9021     Amygdala       BLA         PI+NeuN+
#   4  34ac_scp   Human  Br8331     Amygdala       BLA         PI+NeuN+
#   5  35ac_scp   Human  Br5273     Amygdala       BLA         PI+NeuN+
#   6   Sample1 Macaque    Mac2     Amygdala   Lateral Ventral         

##set up sample data table
sample_info <- data.frame(
    sample_id = paste(tmp$Subject, tmp$Sample, sep = "-"),
    sample_name = tmp$Sample,
    subject = tmp$Subject,
    species = tmp$Species,
    region = tmp$Region,
    subregion = tmp$Subregion,
    dv_axis = tmp$DV_axis
)

stopifnot(all(!duplicated(sample_info$sample_id)))

# add path to cellranger output
sample_info$sample_path<-rep(NA,2)
sample_info$sample_path<- file.path(
    here::here("processed-data", "01_cellranger"),
    "Baboon_ncbi",
    sample_info$sample_name,
    "outs",
    "raw_feature_bc_matrix"
)

# Subset to just baboon samples
sample_info <- sample_info %>%
    filter(species == "Baboon")


## Build basic SCE
message("Read 10x data and create sce - ", Sys.time())
# Read 10x data and create sce - 2023-04-17 13:37:16

sce <- read10xCounts(
    sample_info$sample_path,
    sample_info$sample_name,
    type = "sparse",
    col.names = TRUE
)
message("RDone - ", Sys.time())
# RDone - 2023-04-17 13:38:35

dim(sce)
# [1]    34215 13939066

sum(grepl("^LOC", rownames(sce)))
# [1] 17049

34215-17049
# [1] 17166

# Note that the rownames are the Ensembl gene names - let's switch to the more familiar gene symbols:
# ...but also notice that length(unique(rowData(sce)$Symbol)) != nrow(sce)
#   - That's because some gene symbols are used multiple times.

## Use code from https://github.com/LieberInstitute/Visium_IF_AD/commit/08df3f7e4a3178563d6b4b1861b664b21466b395#diff-10cb35de98e2a3e5f4235cd88f6dabce5469eead2b2db1fd7121126849fcf585L100
## Read in the gene information from the annotation GTF file
gtf <-
    rtracklayer::import(
        "/dcs04/lieber/marmaypag/BLA_crossSpecies_LIBD1070/BLA_crossSpecies/raw-data/refdata/Panubis1_genome_ncbi/genes/genes.gtf"
    )

gtf <- gtf[gtf$type == "gene"]

# find genes ND1, ND2, COX1, COX2
features <- c("ND1", "ND2", "COX1", "COX2")
mt_genes <- gtf[gtf$gene %in% features, ]
mt_genes

# find same genes in sce
mt_genes_sce <- rowData(sce)[rowData(sce)$Symbol %in% features, ]
mt_genes_sce
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id
gtf

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
match_genes
sum(is.na(match_genes))
# [1] 13124

# get only matched genes in SCE
sce <- sce[!is.na(match_genes), ]
dim(sce)
# dim: 21091 13939066 

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "db_xref", "gene", "gene_biotype")]
gtf

## Add the gene info to our SCE object
rowRanges(sce) <- gtf[match_genes]
sce
# class: SingleCellExperiment 
# dim: 21091 13939066 
# metadata(1): Samples
# assays(1): counts
# rownames(21091): TMEM88B ANKRD65 ... ND6 CYTB
# rowData names(3): ID Symbol Type
# colnames(13939066): 1_AAACCCAAGAAACACT-1 1_AAACCCAAGAAACCAT-1 ...
# 7_TTTGTTGTCTTTGGCT-1 7_TTTGTTGTCTTTGTCG-1
# colData names(2): Sample Barcode
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

# add db_xref to gene_ids 
rowData(sce)$gene_id_ncbi <- gsub("GeneID:", "", gtf$db_xref)

# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id_ncbi, rowData(sce)$ID)
rownames(sce) <- rowData(sce)$Symbol.uniq
sce
# class: SingleCellExperiment 
# dim: 21091 13939066 
# metadata(1): Samples
# assays(1): counts
# rownames(21091): TMEM88B ANKRD65 ... ND6 CYTB
# rowData names(5): ID Symbol Type gene_id_ncbi Symbol.uniq
# colnames(13939066): 1_AAACCCAAGAAACACT-1 1_AAACCCAAGAAACCAT-1 ...
# 7_TTTGTTGTCTTTGGCT-1 7_TTTGTTGTCTTTGTCG-1
# colData names(2): Sample Barcode
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):


# Add metadata
sce$key <- paste0(sce$Barcode, "_", sce$Sample)
new_col <- merge(colData(sce), sample_info[, -which(colnames(sample_info) == "sample_path")])
new_col <- new_col[match(sce$key, new_col$key), ]
stopifnot(identical(sce$key, new_col$key))
rownames(new_col) <- colnames(sce)
colData(sce) <- new_col

## Inspect object
sce



# ======== Add MT- to front of NCBI gene names =========
library("singleCellTK")
data("MitoGenes")
MitoGenes

# drop MT- from MitoGenes gene names
MitoGenes.names <- gsub("MT-", "", MitoGenes$human_symbol$human_mito_symbol)

# see if any of the mito genes are in the sce
mt_genes_sce <- rowData(sce)[rowData(sce)$Symbol %in% MitoGenes.names, ]
mt_genes_sce
# ID      Symbol            Type gene_id_ncbi Symbol.uniq
# <character> <character>     <character>  <character> <character>
#     ND1          ND1         ND1 Gene Expression     14444642         ND1
# ND2          ND2         ND2 Gene Expression     14444630         ND2
# ATP8        ATP8        ATP8 Gene Expression     14444633        ATP8
# ATP6        ATP6        ATP6 Gene Expression     14444634        ATP6
# ND3          ND3         ND3 Gene Expression     14444636         ND3
# ND4L        ND4L        ND4L Gene Expression     14444637        ND4L
# ND4          ND4         ND4 Gene Expression     14444638         ND4
# ND5          ND5         ND5 Gene Expression     14444639         ND5
# ND6          ND6         ND6 Gene Expression     14444640         ND6

# add MT- to front of NCBI mito gene names
rowData(sce)$Symbol.uniq <- ifelse(rowData(sce)$Symbol %in% MitoGenes.names, paste0("MT-", rowData(sce)$Symbol), rowData(sce)$Symbol)
# do the same for ID and Symbol
rowData(sce)$ID <- ifelse(rowData(sce)$Symbol %in% MitoGenes.names, paste0("MT-", rowData(sce)$ID), rowData(sce)$ID)
rowData(sce)$Symbol <- ifelse(rowData(sce)$Symbol %in% MitoGenes.names, paste0("MT-", rowData(sce)$Symbol), rowData(sce)$Symbol)


mt_genes_sce <- rowData(sce)[rowData(sce)$ID %in% MitoGenes$human_symbol$human_mito_symbol, ]
mt_genes_sce
# ID      Symbol            Type gene_id_ncbi Symbol.uniq
# <character> <character>     <character>  <character> <character>
#     ND1       MT-ND1      MT-ND1 Gene Expression     14444642      MT-ND1
# ND2       MT-ND2      MT-ND2 Gene Expression     14444630      MT-ND2
# ATP8     MT-ATP8     MT-ATP8 Gene Expression     14444633     MT-ATP8
# ATP6     MT-ATP6     MT-ATP6 Gene Expression     14444634     MT-ATP6
# ND3       MT-ND3      MT-ND3 Gene Expression     14444636      MT-ND3
# ND4L     MT-ND4L     MT-ND4L Gene Expression     14444637     MT-ND4L
# ND4       MT-ND4      MT-ND4 Gene Expression     14444638      MT-ND4
# ND5       MT-ND5      MT-ND5 Gene Expression     14444639      MT-ND5
# ND6       MT-ND6      MT-ND6 Gene Expression     14444640      MT-ND6

rownames(sce) <- rowData(sce)$Symbol.uniq
sce
# class: SingleCellExperiment 
# dim: 21091 13939066 
# metadata(1): Samples
# assays(1): counts
# rownames(21091): TMEM88B ANKRD65 ... MT-ND6 CYTB
# rowData names(5): ID Symbol Type gene_id_ncbi Symbol.uniq
# colnames(13939066): 1_AAACCCAAGAAACACT-1 1_AAACCCAAGAAACCAT-1 ...
# 7_TTTGTTGTCTTTGGCT-1 7_TTTGTTGTCTTTGTCG-1
# colData names(10): Sample Barcode ... subregion dv_axis
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):

if (!dir.exists(here("processed-data", "02_build_sce"))) dir.create(here("processed-data", "02_build_sce"))
save(sce, file = here("processed-data",  "02_build_sce", "sce_baboon_raw_ncbi.rda"))

# get all mito genes. they don't start with MT, so find another way

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
