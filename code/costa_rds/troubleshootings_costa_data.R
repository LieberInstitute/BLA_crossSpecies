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

# load new data from Evan
mat <- readMM(file = "/Users/mtotty2/Documents/R/Costa_NHP/raw-data/samples/CostaData_NonNormalized_Counts.mtx")
length(unique(rownames(mat)))

data <- read.csv(here(samples_dir, "CostaData_CellsAnnos.csv"))
print(data)

unique(data$Batch)





# ========== Read in the CSV file as a data frame ==========
sample.info <- read.delim(here(sampleinfo_dir,"NHP_sampleinfo1.csv"), header = T,sep=',')

baboon_info <- sample.info %>%
  filter(Species == "Baboon")

macaque_info <- sample.info %>%
  filter(Species == "Macaque")



# functions to read cellrange outputs and get gene counts

read_cellranger_outputs <- function(sample) {
  barcode.path <- here(samples_dir, sample, "barcodes.tsv.gz")
  features.path <- here(samples_dir, sample, "features.tsv.gz")
  matrix.path <- here(samples_dir, sample, "matrix.mtx.gz")
  
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  
  return(mat)
}

process_sample_directories <- function(samples) {
  results <- data.frame(Sample = character(), NumGenes = numeric())
  
  for (sample in samples) {
    mat <- read_cellranger_outputs(sample)
    num_genes <- length(unique(rownames(mat)))
    sample_name <- basename(sample)
    results <- rbind(results, data.frame(Sample = sample_name, NumGenes = num_genes))
  }
  
  write.csv(results, here(processed_dir, "sample_gene_counts.csv"), row.names = FALSE)
}


# ==== Baboon data ====
sample_directories <- baboon_info$Sample
process_sample_directories(sample_directories)

# ==== Macaque data ====
sample_directories <- macaque_info$Sample
process_sample_directories(sample_directories)

