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

# Subset to just amygdala macaque samples
sample_info <- sample_info %>%
  filter(species == "Macaque") %>%
  filter(region == "Amygdala")


# Build basic SCE
 message("Read 10x data and create sce - ", Sys.time())
 
 sce <- read10xCounts(
   sample_info$sample_path,
   sample_info$sample_name,
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
gtf
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


# ============= Merging metadata ==============


# Drop sample path before merging metdata
sample_info$sample_path <- NULL

# ===== NOTE =====
# because there are so many samples (~70M droplets), we'll need to use 
# data.tables for more efficient memory as well as data chunking
# to avoid erros from base R's merge() function. 
# ===== NOTE =====

# Convert to data.tables which can more efficiently handle large data
colDataDT <- as.data.table(colData(sce))
sampleInfoDT <- as.data.table(sample_info)

# Ensure the 'Sample' columns are correctly formatted and aligned in both tables
#  (this is the column we're merging on)
colDataDT[, Sample := as.character(Sample)]
sampleInfoDT[, Sample := as.character(Sample)]

# Optionally, if the 'Sample' names are not consistently formatted, you might want to make them consistent
## colDataDT[, Sample := tolower(Sample)]
## sampleInfoDT[, Sample := tolower(Sample)]

# set the number of chunks to be used
chunk_size <- 5000  # This is an arbitrary number; please adjust according to your system's capacity.
number_of_chunks <- ceiling(nrow(colDataDT) / chunk_size)

# Process the data in chunks based on similar columns
mergedResults <- list()
for (i in seq_len(number_of_chunks)) {
  start_row <- ((i - 1) * chunk_size) + 1
  end_row <- min(nrow(colDataDT), i * chunk_size)  # Don't exceed the number of rows
  
  chunk <- colDataDT[start_row:end_row]
  
  # Merge on the 'Sample' column, which exists in both data.tables
  merged_chunk <- merge(chunk, sampleInfoDT, by = "Sample", all.x = TRUE, allow.cartesian = TRUE)
  
  mergedResults[[i]] <- merged_chunk
}

# rbind to concatenate all the chunks back together
# (check here to make sure it looks right)
finalResult <- rbindlist(mergedResults)
# head(finalResult)
# Sample            Barcode                               sample_id subject
# 1: VC_snRNAseq-7_LateralVentralAP1_-1 AAACCCAAGAAACCCA-1 Mac2-VC_snRNAseq-7_LateralVentralAP1_-1    Mac2
# 2: VC_snRNAseq-7_LateralVentralAP1_-1 AAACCCAAGAAACCCG-1 Mac2-VC_snRNAseq-7_LateralVentralAP1_-1    Mac2
# 3: VC_snRNAseq-7_LateralVentralAP1_-1 AAACCCAAGAAACTAC-1 Mac2-VC_snRNAseq-7_LateralVentralAP1_-1    Mac2
# 4: VC_snRNAseq-7_LateralVentralAP1_-1 AAACCCAAGAAACTCA-1 Mac2-VC_snRNAseq-7_LateralVentralAP1_-1    Mac2
# 5: VC_snRNAseq-7_LateralVentralAP1_-1 AAACCCAAGAAAGTCT-1 Mac2-VC_snRNAseq-7_LateralVentralAP1_-1    Mac2
# 6: VC_snRNAseq-7_LateralVentralAP1_-1 AAACCCAAGAAATGAG-1 Mac2-VC_snRNAseq-7_LateralVentralAP1_-1    Mac2
# species   region subregion dv_axis
# 1: Macaque Amygdala   Lateral Ventral
# 2: Macaque Amygdala   Lateral Ventral
# 3: Macaque Amygdala   Lateral Ventral
# 4: Macaque Amygdala   Lateral Ventral
# 5: Macaque Amygdala   Lateral Ventral
# 6: Macaque Amygdala   Lateral Ventral

# Convert 'finalResult' from a data.table/data.frame to a DataFrame
finalResultDF <- DataFrame(finalResult)

# Reordering if needed (assuming 'cell_id' is your unique cell identifier column)
finalResultDF <- finalResultDF[match(colData(sce)$Barcode, finalResultDF$Barcode),]

# Completely replace the existing colData
colData(sce) <- finalResultDF

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