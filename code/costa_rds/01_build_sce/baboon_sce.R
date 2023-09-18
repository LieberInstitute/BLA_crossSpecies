library("SingleCellExperiment")
library("DropletUtils")
library("here")
library("rtracklayer")
library("sessioninfo")
library("dplyr")


# directories
raw_dir = here("raw-data")
samples_dir = here(raw_dir, "samples")
sampleinfo_dir = here(raw_dir, "sampleinfo")
processed_dir = here("processed-data")

# ========== Read in the CSV file as a data frame ==========
tmp <- read.delim(here(sampleinfo_dir,"NHP_sampleinfo.csv"), header = T,sep=',')

baboon_info <- tmp %>%
  filter(Species == "Baboon")

# View the subsetted data
baboon_info
# 1         Bbn2_Basal1   Basal            Baboon2   Male  Baboon    NA         Sample1
# 2         Bbn2_Basal2   Basal            Baboon2   Male  Baboon    NA         Sample2
# 3       Bbn2_Lateral1 Lateral            Baboon2   Male  Baboon    NA         Sample3
# 4         Bbn5_Basal1   Basal            Baboon5 Female  Baboon    NA         Sample4
# 5    Bbn5_Basal1_Plus   Basal            Baboon5 Female  Baboon    NA         Sample6
# 6       Bbn5_Lateral1 Lateral            Baboon5 Female  Baboon    NA         Sample7
# 7  Bbn5_Lateral1_Plus Lateral            Baboon5 Female  Baboon    NA         Sample5
# 8           An3_BA_L1   Basal                    Female  Baboon    NA         Sample8
# 9           An3_LA_L1 Lateral                    Female  Baboon    NA         Sample9
# 10          An3_LA_L2 Lateral                    Female  Baboon    NA        Sample10

# ========== Set up sample data frame ==========
sample_info <- data.frame(
  #sample_id = paste(tmp$Subject, tmp$Sample, sep = "-"),
  Sample = baboon_info$Sample,
  Region = baboon_info$Region,
  Sex = baboon_info$Sex,
  Anatomoical = baboon_info$Anatomical,
  Batch = baboon_info$Batch
)

head(sample_info)
# Sample  Region    Sex Anatomoical Batch
# 1      Baboon2_Basal1   Basal   Male                NA
# 2      Baboon2_Basal2   Basal   Male                NA
# 3    Baboon2_Lateral1 Lateral   Male                NA
# 4      Baboon5_Basal1   Basal Female                NA
# 5 Baboon5_Basal1_Plus   Basal Female                NA
# 6    Baboon5_Lateral1 Lateral Female                NA

# check for duplicates
stopifnot(all(!duplicated(sample_info$Sample)))

# add path to cellranger output
sample_info$sample_path
sample_info$sample_path<- file.path(
  samples_dir,
  sample_info$Sample
)

sample_info$sample_path
# [1] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/Bbn2_Basal1"     
# [2] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/Bbn2_Basal2"     
# [3] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/Bbn2_Lateral1"   
# [4] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/Bbn5_Basal1"     
# [5] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/Bbn5_Basal1_Plus"
# [6] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/Bbn5_Lateral1" 
# [6] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/Bbn5_Lateral1"     
# [7] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/Bbn5_Lateral1_Plus"
# [8] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/An3_BA_L1"         
# [9] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/An3_LA_L1"         
# [10] "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/An3_LA_L2" 

# ============= Build basic SCE ==============
# I initially tried to run this with both Macacue and Baboon samples and received
# the error "gene information differs between runs".

# After running seperately, it runs fine.
message("Read 10x data and create sce - ", Sys.time())
# Read 10x data and create sce - 2023-08-03 07:42:28


sce <- read10xCounts(
  sample_info$sample_path,
  sample.names = sample_info$Sample,
  type = "sparse",
  col.names = TRUE
)
message("RDone - ", Sys.time())
# RDone - 2023-08-03 07:50:18


# ===== PROBLEM ======
# it seems like a lot of the "genes" (~16k) are these LOC names, which do not map
# to actual genes in the baboon. Could be due to inappropriate genome matching?
> head(grep("LOC", rownames(sce), value = TRUE))
# [1] "LOC101011345" "LOC103884675" "LOC116273776" "LOC116273781" "LOC103881263" "LOC103881261"



# ============= Unique-ify gene names ============= 
# Note that the rownames are the Ensembl gene names - let's switch to the more familiar gene symbols:
# ...but also notice that length(unique(rowData(sce)$Symbol)) != nrow(sce)
#   - That's because some gene symbols are used multiple times.

## Use files from https://ftp.ensembl.org/pub/release-110/gtf/macaca_mulatta/

# gtf <-
#   rtracklayer::import(
#     "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
#   )

gtf <-
  rtracklayer::import(
    "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/genes/Panubis1.0.110.gtf"
  )

gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_name

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_name)
matched_genes <- match_genes[!is.na(match_genes)]

# The genes were not a match.
# I used the human genes - so obviously they wouldn't match. 
# Will get macacue genes and update. 


## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_biotype")]

## Add the gene info to our SPE object
rowRanges(sce) <- gtf[match_genes]

# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq


# ========== Add metadata ===========
# backup.sce <- sce

sce$sample_id <- sce$Sample

sce$key <- paste0(sce$Barcode, "_", sce$Sample)
new_col <- merge(colData(sce), sample_info)
new_col <- new_col[match(sce$key, new_col$key), ]
stopifnot(identical(sce$key, new_col$key))
rownames(new_col) <- colnames(sce)
colData(sce) <- new_col

## Inspect object
sce
# class: SingleCellExperiment 
# dim: 22353 152472 
# metadata(1): Samples
# assays(1): counts
# rownames(22353): U6_ENSMMUG00000036181 ZNF692 ... UBE2M KIR3DH
# rowData names(7): source type ... gene_biotype Symbol.uniq
# colnames(152472): 1_AAACCCAAGAACGCGT-1 1_AAACCCAAGCACCAGA-1 ... 41_TTTGGTTGTGGCTAGA-1 41_TTTGTTGCACATGGTT-1
# colData names(8): Sample Barcode ... Batch sample_path
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):

unique(sce$Batch)
# [1] 1 2 3 4 5

unique(sce$Sample)
# [1] "BA_1"                 "BA_2"                 "LA_1"                 "LA_2"                 "AB_1A_1"             
# [6] "AB_1A_2"              "AB_1B_1"              "AB_1B_2"              "LateralVentralAP1"    "LateralVentralAP2"   
# [11] "BasalDorsalAP1_3A"    "BasalDorsalAP1_3B"    "BasalVentralAP1"      "BasalVentralAP2"      "AccBasalAP1AP2"      
# [16] "LateralDorsalAP1AP2"  "LV_13"                "LV3_13"               "LD_13"                "BV2_13"              
# [21] "BV1_13"               "BD_13"                "AB_13"                "Bd_Bv_13"             "LV2_14"              
# [26] "LV1_14"               "LD_14"                "L_Comb_14"            "BV_14"                "BD_14"               
# [31] "B_Comb_14"            "AB_14"                "CN_S7"                "CN_S8"                "CN_S1"               
# [36] "AccessoryBasalAP2AP3" "BasalAP1_S2"          "BasalAP2_S5"          "Lateral_AP1_S4"       "Lateral_AP2_S6"      
# [41] "Central_Nucleus_AP3" 

unique(sce$Anatomoical)
# [1] "Ventral"        ""               "Dorsal"         "Dorsal_Ventral"

## Size in Gb
lobstr::obj_size(sce) / 1024^3
# 4.44 B

if (!dir.exists(here("processed-data", "01_build_sce"))) dir.create(here("processed-data", "01_build_sce"))
save(sce, file = here("processed-data", "01_build_sce", "sce_maqacue.rda"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()