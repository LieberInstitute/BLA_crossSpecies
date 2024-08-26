library(here)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(dplyr)


# directories
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "08_species_comparisons")


sce.excit <- readRDS(here("processed-data", "07_annotation_and_characterization", "sce_excit_final_subclusters_annotated.rds"))
sce.inhib <- readRDS(here("processed-data", "07_annotation_and_characterization", "sce_inhib_final_subclusters_annotated.rds"))




# ==== add putative subregion labels ====

# combine putative LA, BA, and aBA clusters
LA_clusters <- c("GULP1_TRHDE", "ZBTB20_SLC4A4")
BA_clusters <- c("PEX5L_MYRIP", "MEIS2_COL25A1")
aBA_clusters <- c("ESR1_ADRA1A", "GRIK3_TNS3")

sce.excit$subregion_celltype <- "Other"
sce.excit$subregion_celltype[sce.excit$fine_celltype %in% LA_clusters] <- "LA"
sce.excit$subregion_celltype[sce.excit$fine_celltype %in% BA_clusters] <- "BA"
sce.excit$subregion_celltype[sce.excit$fine_celltype %in% aBA_clusters] <- "aBA"
unique(sce.excit$subregion_celltype)
# [1] "aBA"   "BA"    "Other" "LA"   


# genes <- c(
#     "GRM1", "GRM5", "GRM2", "GRM3", "GRM4", "GRM6", "GRM7", "GRM8",
#     "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5",
#     "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B",
#     "HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR2A", "HTR2C", "HTR4", "HTR5A", "HTR6", "HTR7",
#     "DRD1", "DRD2", "DRD3", "DRD4", "DRD5",
#     "CHRNA3", "CHRNA5", "CHRNA6", "CHRNB3", "CHRNB4",
#     "ADRA1A", "ADRA1B", "ADRA2A", "ADRA2C", "ADRB1", "ADRB2", "AVPR1A", "AVPR1B", "AVPR2", "CNR1", "CNR2", 
#     "DRD2", "DRD3", "F2RL1", "F2RL2", "F2RL3", "GHRHR", "HCRTR2", "KISS1R", "MC2R", "MC3R", "MC4R", "MC5R", 
#     "NPBWR1", "NPBWR2", "NPY1R", "NPY2R", "NPY4R", "NPY5R", "OPRD1", "OPRK1", "OPRM1", "OXTR", "PRLHR", 
#     "PTGER1", "PTGER2", "PTGER3", "PTGER4", "PTH1R", "PTH2R", "SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5", 
#     "TACR1", "TACR2", "TACR3", "TRHR", "V1AR", "V1BR"
# )


serotonin_genes <- c("HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR2A", "HTR2C", "HTR4", "HTR5A", "HTR6", "HTR7")

dopamine_genes <- c("DRD1", "DRD2", "DRD3", "DRD4", "DRD5")

norepinephrine_genes <- c("ADRA1A", "ADRA1B", "ADRA2A", "ADRA2C", "ADRB1", "ADRB2")

cannabinoid_genes <- c("CNR1", "CNR2")

cholinergic_genes <- c("CHRNA3", "CHRNA5", "CHRNA6", "CHRNB3", "CHRNB4")



# drop any not in sce.excit
genes <- genes[genes %in% rownames(sce.excit)]

# drop any duplicate
genes <- unique(genes)

# excitatory
png(here(plot_dir, "dotplot_excitatory_receptrs_and_such.png"), width = 20, height = 5, units = "in", res = 300)
plotDots(sce.excit, features=serotonin_genes, group="fine_celltype", block="species", center=TRUE, scale=TRUE) +
    coord_flip() +
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    # x label horizontal
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
     facet_wrap("subregion_celltype"~., scales = "free_y") 
dev.off()


# inhibitory
png(here(plot_dir, "dotplot_inhibitory_receptrs_and_such.png"), width = 20, height = 5, units = "in", res = 300)
plotDots(sce.inhib, features=genes, group="fine_celltype", block="species", center=TRUE, scale=TRUE) +
    coord_flip() +
    scale_color_gradient2(low = "blue", mid="white", high = "red") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(color="centered.scaled") 
dev.off()




# Aggregate mean expression across cells by fine_celltype
agg_data <-scuttle::aggregateAcrossCells(
    sce.excit,
    ids = sce.excit$fine_celltype, # Aggregate by fine_celltype
    statistics = "mean",
    assay.type="logcounts",
    subset.row=serotonin_genes
)



# Create a data frame for gene categories
gene_categories <- data.frame(
    gene = c(serotonin_genes, dopamine_genes, norepinephrine_genes, cannabinoid_genes, cholinergic_genes),
    category = c(rep("Serotonin", length(serotonin_genes)),
                 rep("Dopamine", length(dopamine_genes)),
                 rep("Norepinephrine", length(norepinephrine_genes)),
                 rep("Cannabinoid", length(cannabinoid_genes)),
                 rep("Cholinergic", length(cholinergic_genes)))
)

# Add the expression data from the sce.excit object (assuming it's in a tidy format)





