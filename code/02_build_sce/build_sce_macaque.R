library("SingleCellExperiment")
library("DropletUtils")
library("here")
library("rtracklayer")
library("sessioninfo")
library("dplyr")
library("data.table")

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
  Sample = tmp$Sample,
  subject = tmp$Subject,
  species = tmp$Species,
  region = tmp$Region,
  subregion = tmp$Subregion,
  dv_axis = tmp$DV_axis
)
head(sample_info)

stopifnot(all(!duplicated(sample_info$sample_id)))

# add path to cellranger output
sample_info$sample_path<-rep(NA,2)
sample_info$sample_path<- file.path(
  here::here("processed-data", "01_cellranger"),
  sample_info$species,
  sample_info$Sample,
  "outs",
  "raw_feature_bc_matrix"
)
head(sample_info)

# Subset to just amygdala macaque samples
sample_info <- sample_info %>%
  filter(species == "Macaque") %>%
  filter(region == "Amygdala")


# Build basic SCE
 message("Read 10x data and create sce - ", Sys.time())
 
 sce <- read10xCounts(
   sample_info$sample_path,
   sample_info$Sample,
   type = "sparse",
   col.names = TRUE
 )
 message("RDone - ", Sys.time())

# Note that the rownames are the Ensembl gene names - let's switch to the more familiar gene symbols:
# ...but also notice that length(unique(rowData(sce)$Symbol)) != nrow(sce)
#   - That's because some gene symbols are used multiple times.



#load(here('processed-data', '01_build_sce', 'sce_macaque_raw.rda'))
sce

## Use code from https://github.com/LieberInstitute/Visium_IF_AD/commit/08df3f7e4a3178563d6b4b1861b664b21466b395#diff-10cb35de98e2a3e5f4235cd88f6dabce5469eead2b2db1fd7121126849fcf585L100
## Read in the gene information from the annotation GTF file
gtf <-
  rtracklayer::import(
    here('raw-data', 'refdata' ,'Mmul_10_genome', 'genes', 'genes.gtf')
  )

gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_biotype")]

## Add the gene info to our SPE object
rowRanges(sce) <- gtf[match_genes]


# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq
sce

colnames(colData(sce))

# ============= Merging metadata ==============



# Create the key column for sce
sce$key <- paste0(sce$Barcode, "_", sce$Sample)

# Convert your data to data.table objects
sce_colData <- as.data.table(colData(sce))
sample_info_dt <- as.data.table(sample_info)

# Remove the 'sample_path' column from sample_info_dt
sample_info_dt <- sample_info_dt[, -which(names(sample_info_dt) == "sample_path"), with = FALSE]

# Identify the common columns to merge on
common_cols <- intersect(names(sce_colData), names(sample_info_dt))

# Perform the merge using the [.data.table method
new_col <- sample_info_dt[sce_colData, on = common_cols]

new_col <- new_col[match(sce$key, sce_colData$key), ]

# Check that the keys match
stopifnot(identical(sce$key, new_col$key))

# Convert new_col to a DataFrame as required by SingleCellExperiment
new_col_df <- DataFrame(new_col)

# Update the colData
colData(sce) <- new_col_df

# Optionally, update the row names if needed
rownames(sce) <- rownames(new_col_df)


## Inspect object
sce
# class: SingleCellExperiment 
# dim: 21369 70879726 
# metadata(1): Samples
# assays(1): counts
# rownames(21369): ENSMMUG00000023296 ZNF692 ... ND6 CYTB
# rowData names(7): source type ... gene_biotype Symbol.uniq
# colnames: NULL
# colData names(8): Sample Barcode ... subregion dv_axis
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):

if (!dir.exists(here("processed-data", "02_build_sce"))) dir.create(here("processed-data", "02_build_sce"))
save(sce, file = here("processed-data",  "02_build_sce", "sce_macaque_raw.rda"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()