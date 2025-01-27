library(dplyr)
library(here)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(tidyr)
library(scater)
library(dreamlet)
library(variancePartition)
library(edgeR)


# directories
processed_dir = here("processed-data", "07_annotation_and_characterization")
plot_dir = here("plots", "08_species_comparisons", "03_DE_analysis","dream")

sce.excit <- readRDS(here(processed_dir, "sce_excit_final_subclusters_annotated.rds"))


# ======== Pseudobulk analysis ========
# for DE analysis, it makes the most sense to subset to anatomical samples
# and perform pseudobulk and DE analysis. This avoids the statistical "double dipping"
# of using markers genes based on clusters... which are based on DEGs.

# subset SCE to just Later and Basal
sce.subset <- sce.excit[, sce.excit$Subregion %in% c("Lateral", "Basal")]
dim(sce.subset)
# [1] 13874 53306



# aggregate across cells
pseudo <- aggregateAcrossCells(sce.subset, ids = colData(sce.subset)[,c("Subregion", "species", "Sample","Subject")])
dim(pseudo)
#[1] 13874    32




dge <- DGEList(cpm(pseudo))
dge <- calcNormFactors(dge)

metadata<- colData(pseudo)
# ====== DE using Dreamlet =======

# The variable to be tested must be a fixed effect
form <- ~ Subregion  + (1 | Subject) + (1 | species)

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(dge, form, metadata)

# Fit the dream model on each gene
# For the hypothesis testing, by default,
# dream() uses the KR method for <= 20 samples,
# otherwise it uses the Satterthwaite approximation
fitmm <- dream(vobjDream, form, metadata)
fitmm <- eBayes(fitmm)