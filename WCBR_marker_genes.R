library(dplyr)
library("Matrix")
library("here")
library("ggspavis")
library("scater")
library("pheatmap")
library("spatialLIBD")
library("patchwork")
library("scran")
library("Seurat")
library("SingleCellExperiment")
library(pheatmap)
library(MetaNeighbor)
library(tidyverse)
library(ggplot2)
library(ghibli)
library(RColorBrewer)
library(ComplexHeatmap)

# directories
processed_dir = here("processed-data", "07_annotation")
plot_dir = here("plots", "07_annotation", "WCBR")

sce.inhib <- readRDS(here(processed_dir, "sce.inhib.final.rds"))
sce.excit <- readRDS(here(processed_dir, "sce.excit.integrated.annotated.rds"))



# ======== Inhibitory markers ========


# ===== marker genes =====
features <- c("CARTPT", "CCK", "LAMP5",
              "PENK", "PPP1R1B",
               "PVALB", "RELN", "SST", "NOS1",
              "TSHZ1", "VIP", "ZFHX3"
)

# drop features not in sce.inhib
features <- features[features %in% rownames(sce.inhib)]
features

pdf(here(plot_dir, "Dotplot_inhib_marker_genes.pdf"), width=7.5, height=7.5)
plotDots(sce.inhib, 
         features = features, 
         group = "ident", 
         center=TRUE, 
         scale=TRUE)
dev.off()

pdf(here(plot_dir, "Heatmap_inhib_marker_genes_top.pdf"), width=5, height=5)
scater::plotGroupedHeatmap(sce.inhib, 
            features = features,
            group = "ident", 
            block="species",
            center=TRUE, 
            scale=TRUE,
            labels_col=unique(sce.inhib$fine_type),
            cluster_rows=FALSE,
            cluster_cols=FALSE
            )
dev.off()


# trying to make this better by plotting along the diagonal using sheatmap
library(slanter)

sheatmap(sce.inhib, 
         features = features,
         group = "ident", 
         block="species",
         center=TRUE, 
         scale=TRUE,
         labels_col=unique(sce.inhib$fine_type),
         cluster_rows=FALSE,
         cluster_cols=FALSE
)



out <- aggregateAcrossCells(sce.inhib,
                            ids=sce.inhib$fine_type, 
                            subset.row = features,
                            use.assay.type = "logcounts"
                            )

df <- data.frame(logcounts(out))

pdf(here(plot_dir, "Heatmap_inhib_marker_genes.pdf"), width=5, height=5)
sheatmap(data.matrix(df),
         #center=TRUE,
         scale="row",
         order_cols=FALSE,
         order_rows=FALSE,
         treeheight_row = 0, 
         treeheight_col = 0
         )
dev.off()
# ======== Excitatory markers ========

features <- c("COL25A1",
              "GULP1",
              "ZNF804B",
              "TENM3",
              "CHRM3",
              "SAMD5", 
              "GRM8"
              #"GPC5"
)

# drop features not in sce.inhib
features <- features[features %in% rownames(sce.excit)]
features

pdf(here(plot_dir,"Dotplot_excit_marker_genes_top.pdf"), width=7.5, height=7.5)
plotDots(sce.excit, 
         features = features, 
         group = "ident", 
         center=TRUE, 
         scale=TRUE)
dev.off()




out <- aggregateAcrossCells(sce.excit,
                            ids=sce.excit$fine_type, 
                            subset.row = features,
                            use.assay.type = "logcounts"
)

df <- data.frame(logcounts(out))

pdf(here(plot_dir, "Heatmap_excit_marker_genes.pdf"), width=5, height=5)
sheatmap(data.matrix(df),
         center=TRUE,
         scale="column",
         #order_cols=FALSE,
         #order_rows=FALSE,
         treeheight_row = 0, 
         treeheight_col = 0
)
dev.off()