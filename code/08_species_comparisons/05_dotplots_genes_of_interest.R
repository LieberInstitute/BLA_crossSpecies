library(here)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(scater)


# directories
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "08_species_comparisons")


sce.excit <- readRDS(here("processed-data", "sce_excit_final_subclusters_annotated.rds"))
sce.inhib <- readRDS(here("processed-data", "sce_inhib_final_subclusters_annotated.rds"))




genes <- c(
    "GRM1", "GRM5", "GRM2", "GRM3", "GRM4", "GRM6", "GRM7", "GRM8",
    "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5",
    "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B",
    "HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR2A", "HTR2C", "HTR4", "HTR5A", "HTR6", "HTR7",
    "DRD1", "DRD2", "DRD3", "DRD4", "DRD5",
    "CHRNA3", "CHRNA5", "CHRNA6", "CHRNB3", "CHRNB4",
    "ADRA1A", "ADRA1B", "ADRA2A", "ADRA2C", "ADRB1", "ADRB2", "AVPR1A", "AVPR1B", "AVPR2", "CNR1", "CNR2", 
    "DRD2", "DRD3", "F2RL1", "F2RL2", "F2RL3", "GHRHR", "HCRTR2", "KISS1R", "MC2R", "MC3R", "MC4R", "MC5R", 
    "NPBWR1", "NPBWR2", "NPY1R", "NPY2R", "NPY4R", "NPY5R", "OPRD1", "OPRK1", "OPRM1", "OXTR", "PRLHR", 
    "PTGER1", "PTGER2", "PTGER3", "PTGER4", "PTH1R", "PTH2R", "SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5", 
    "TACR1", "TACR2", "TACR3", "TRHR", "V1AR", "V1BR"
)

# drop any not in sce.excit
genes <- genes[genes %in% rownames(sce.excit)]

# drop any duplicate
genes <- unique(genes)

# excitatory
png(here(plot_dir, "dotplot_excitatory_receptrs_and_such.png"), width = 20, height = 5, units = "in", res = 300)
plotDots(sce.excit, features=genes, group="fine_celltype", block="species", center=TRUE, scale=TRUE) +
    coord_flip() +
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    # x label horizontal
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()


# inhibitory
png(here(plot_dir, "dotplot_inhibitory_receptrs_and_such.png"), width = 20, height = 5, units = "in", res = 300)
plotDots(sce.inhib, features=genes, group="fine_celltype", block="species", center=TRUE, scale=TRUE) +
    coord_flip() +
    scale_color_gradient2(low = "blue", mid="white", high = "red") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(color="centered.scaled")
dev.off()








