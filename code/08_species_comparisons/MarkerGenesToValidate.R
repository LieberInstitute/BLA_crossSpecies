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


# subset sce.excit to humans
sce.human <- sce.excit[, sce.excit$species == "human"]
sce.nhp <- sce.excit[, sce.excit$species != "human"]

# Lateral markers
features <- c("ZBTB20", "ZNF536","GULP1", "HS6ST3", "TSHZ2", "CDH6", "NCALD", "GABRB2", "KCTD8", "AVP")
png(here(plot_dir, "lateral_markers_human.png"), width=10, height=10, units="in", res=300)
plotExpression(sce.human, features=features, x="fine_type", color="fine_type", ncol=1) +
    ggtitle("Human: Lateral markers")
dev.off()

png(here(plot_dir, "lateral_markers_nhp.png"), width=10, height=10, units="in", res=300)
plotExpression(sce.nhp, features=features, x="fine_type", color="fine_type", ncol=1) +
    ggtitle("NHP: Lateral markers")
dev.off()

# Basal markers
features <- c("SLIT2", "COL25A1","DPP10", "MYRIP", "PEX5L", "LAMP5", "RGS6")
png(here(plot_dir, "basal_markers_human.png"), width=10, height=10, units="in", res=300)
plotExpression(sce.human, features=features, x="fine_type", color="fine_type", ncol=1) +
    ggtitle("Human: Lateral markers")
dev.off()

png(here(plot_dir, "basal_markers_nhp.png"), width=10, height=10, units="in", res=300)
plotExpression(sce.nhp, features=features, x="fine_type", color="fine_type", ncol=1) +
    ggtitle("NHP: Lateral markers")
dev.off()

