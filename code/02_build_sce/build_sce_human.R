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
  sample_info$species,
  sample_info$sample_name,
  "outs",
  "raw_feature_bc_matrix"
)

# Subset to just human samples
sample_info <- sample_info %>%
  filter(species == "Human")


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

# Note that the rownames are the Ensembl gene names - let's switch to the more familiar gene symbols:
# ...but also notice that length(unique(rowData(sce)$Symbol)) != nrow(sce)
#   - That's because some gene symbols are used multiple times.

## Use code from https://github.com/LieberInstitute/Visium_IF_AD/commit/08df3f7e4a3178563d6b4b1861b664b21466b395#diff-10cb35de98e2a3e5f4235cd88f6dabce5469eead2b2db1fd7121126849fcf585L100
## Read in the gene information from the annotation GTF file
gtf <-
    rtracklayer::import(
        "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our SPE object
rowRanges(sce) <- gtf[match_genes]

# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq


# Add metadata
sce$key <- paste0(sce$Barcode, "_", sce$Sample)
new_col <- merge(colData(sce), sample_info[, -which(colnames(sample_info) == "sample_path")])
new_col <- new_col[match(sce$key, new_col$key), ]
stopifnot(identical(sce$key, new_col$key))
rownames(new_col) <- colnames(sce)
colData(sce) <- new_col

## Inspect object
sce



if (!dir.exists(here("processed-data", "02_build_sce"))) dir.create(here("processed-data", "02_build_sce"))
save(sce, file = here("processed-data", "02_build_sce", "sce_human_raw.rda"))


## Size in Gb
#lobstr::obj_size(sce) / 1024^3
# 1.27

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

