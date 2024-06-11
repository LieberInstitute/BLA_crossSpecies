library(dplyr)
library(here)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(tidyr)
library(scater)


# directories
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "08_species_comparisons", "03_DE_analysis")

sce.inhib <- readRDS(here(processed_dir, "sce.inhib.final.rds"))
sce.excit <- readRDS(here(processed_dir, "sce.excit.integrated.annotated.rds"))





# ======== Pseudobulk analysis ========
# for DE analysis, it makes the most sense to subset to anatomical samples
# and perform pseudobulk and DE analysis. This avoids the statistical "double dipping"
# of using markers genes based on clusters... which are based on DEGs.

# subset SCE to just Later and Basal
sce.subset <- sce.excit[, sce.excit$Subregion %in% c("Lateral", "Basal")]
dim(sce.subset)
# [1] 13874 53306

# aggregate across cells
pseudo <- aggregateAcrossCells(sce.subset, ids = colData(sce.subset)[,c("Subregion", "species","DV_axis", "Sample")])
dim(pseudo)
#[1] 13874    32

unique(pseudo$species)
# [1] "baboon"  "macaque"

# filter low expressing genes
keep <- edgeR::filterByExpr(pseudo, group=pseudo$Subregion)
pseudo <- pseudo[keep,]
dim(pseudo)



library(DeconvoBuddies)

ratio_markers <- get_mean_ratio2(
    sce.subset,
    cellType_col = "Subregion",
    assay_name = "logcounts",
    add_symbol = TRUE
)


top_ratio_markers <- ratio_markers %>%
    group_by(cellType.target) %>%
    top_n(n = 100, wt = ratio) %>%
    filter(mean.target > 1)

write.csv(top_ratio_markers, here(processed_dir, "crossSpecies_top_100_ratio_markers_LAvBA.csv"), row.names = FALSE)

markers_1vAll <- findMarkers_1vAll(sce.subset,
                                     assay_name = "logcounts",
                                     cellType_col = "Subregion",
                                     add_symbol = FALSE,
                                     mod = "~species",
                                     verbose = TRUE)

top_markers <- markers_1vAll %>%
    group_by(cellType.target) %>%
    top_n(n = 100, wt = anno_logFC)

write.csv(top_markers, here(processed_dir, "crossSpecies_top_100_markers_LAvBA.csv"), row.names = FALSE)
