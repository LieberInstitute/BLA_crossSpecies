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

# directories
raw_dir = here("raw-data")
samples_dir = here(raw_dir, "samples")
sampleinfo_dir = here(raw_dir, "sampleinfo")
processed_dir = here("processed-data")
plot_dir = here("plots", "05_species_comparisons")

Amy.700 <- readRDS(here(processed_dir, "00_costa_rds", "amygdalaCNScaled700.rds"))
Amy.700
sce.costa <- as.SingleCellExperiment(Amy.700)
# An object of class Seurat 
# 21353 features across 120794 samples within 2 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

# Bab.200 <-readRDS(here(processed_dir, "00_costa_rds", "AllBaboon_Amygdala200_04.rds"))
# Bab.200
# sce.bab <- as.SingleCellExperiment(Bab.200)

human.yu <- readRDS(here(processed_dir, "00_costa_rds", "GSE195445_Human_obj.rda"))
human.yu
sce.yu <- as.SingleCellExperiment(human.yu)

# remove "Human" form the beginning of $ident
levels(sce.yu$ident) <- gsub("Human_", "", levels(sce.yu$ident))

load(here(processed_dir, "00_costa_rds", "sce_annotated.rda"))
sce.totty <- sce



# ====== Subset inhibitory celltypes ======

sce.totty <- sce.totty[, sce.totty$celltype == "Inhib"]
sce.yu <- sce.yu[, sce.yu$celltype == "InN"]


# ====== MetaNeighbor ======

common_genes <- intersect(rownames(sce.totty), rownames(sce.yu))
sce.totty <- sce.totty[common_genes,]
sce.yu <- sce.yu[common_genes, ]

new_colData = data.frame(
  study_id = rep(c('Yu', 'Totty'), c(ncol(sce.yu), ncol(sce.totty))),
  cell_type = c(as.character(colData(sce.yu)$ident), colData(sce.totty)$annotation)
)
sce.bla <- SingleCellExperiment(
  Matrix(cbind(assay(sce.yu, 1), assay(sce.totty, 1)), sparse = TRUE),
  colData = new_colData
)
dim(sce.bla)
rm(sce.totty); rm(sce.yu)

# ====== run it ======

var_genes = variableGenes(dat = sce.bla, exp_labels = sce.bla$study_id)
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = sce.bla,
                             study_id = sce.bla$study_id,
                             cell_type = sce.bla$cell_type,
                             fast_version = TRUE)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
gplots::heatmap.2(celltype_NV,
                  margins=c(8,8),
                  keysize=1,
                  key.xlab="AUROC",
                  key.title="NULL",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 0.7,
                  cexCol = 0.7)

# ========== Identify top hits ==========

top_hits = topHits(cell_NV = celltype_NV,
                   dat = sce.bla,
                   study_id = sce.bla$study_id,
                   cell_type = sce.bla$cell_type,
                   threshold = 0.9)
top_hits
# Study_ID|Celltype_1    Study_ID|Celltype_2 Mean_AUROC         Match_type
# 1     Yu|PVALB ADAMTS5         Totty|Inh_5       1.00 Reciprocal_top_hit
# 2     Yu|LAMP5 COL14A1   Totty|Inh_LAMP5_2       1.00 Reciprocal_top_hit
# 3        Yu|HTR3A DRD2         Totty|Inh_4       1.00 Reciprocal_top_hit
# 4        Yu|LAMP5 NDNF   Totty|Inh_LAMP5_1       1.00 Reciprocal_top_hit
# 5      Yu|TSHZ1 SEMA3C     Totty|Inh_ITC_2       0.99 Reciprocal_top_hit
# 6      Yu|TSHZ1 CALCRL     Totty|Inh_ITC_1       0.98 Reciprocal_top_hit
# 7          Yu|VIP NDNF         Totty|Inh_2       0.98 Reciprocal_top_hit
# 8          Yu|SST EPYC       Totty|Inh_SST       0.98 Reciprocal_top_hit
# 9           Yu|SST HGF     Totty|Inh_PVALB       0.98 Reciprocal_top_hit
# 10         Yu|VIP NDNF         Totty|Inh_1       0.98          Above_0.9
# 11       Yu|VIP ABI3BP     Totty|Inh_VIP_1       0.98 Reciprocal_top_hit
# 12       Yu|VIP ABI3BP     Totty|Inh_VIP_2       0.95          Above_0.9
# 13       Yu|VIP ABI3BP         Totty|Inh_3       0.92          Above_0.9


# backup
sce.backup <- sce.bla



# ====== Run AUROC ======
sce.bla <- sce.backup

# Rename cell types in sce.bla based on the above table
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_5"] <- "PVALB ADAMTS5"
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_LAMP5_2"] <- "LAMP5 COL14A1"
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_4"] <- "HTR3A DRD2"
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_LAMP5_1"] <- "LAMP5 NDNF"
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_ITC_2"] <- "TSHZ1 SEMA3C"
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_ITC_1"] <- "TSHZ1 CALCRL"
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_2"] <- "VIP NDNF"
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_SST"] <- "SST EPYC"
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_PVALB"] <- "SST HGF"
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_1"] <- "VIP NDNF"
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_VIP_1"] <- "VIP ABI3BP"
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_VIP_2"] <- "VIP ABI3BP"
colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == "Inh_3"] <- "VIP ABI3BP"


celltypes_of_interest = sce.bla[, sce.bla$cell_type %in% c('TSHZ1 SEMA3C', 'PVALB ADAMTS5', 'VIP ABI3BP',
                                                           'LAMP5 COL14A1', 'HTR3A DRD2', 'LAMP5 NDNF', 
                                                           'TSHZ1 CALCRL', 'VIP NDNF', 'SST EPYC', 
                                                           'SST HGF', 'VIP ABI3BP')]
celltype_matrix = model.matrix(~celltypes_of_interest$cell_type - 1)
colnames(celltype_matrix) = levels(as.factor(celltypes_of_interest$cell_type))

#data(GOhuman)
pdf(here(plot_dir,"Metaneighbor_AUROC_interneurons.pdf"))
AUROC_scores = MetaNeighbor(dat = celltypes_of_interest ,
                            experiment_labels = as.numeric(factor(celltypes_of_interest$study_id)),
                            celltype_labels = celltype_matrix,
                            genesets = GOhuman,
                            bplot = TRUE,
                            fast_version = TRUE)
dev.off()

# plotting using ggplot
long_data <- data.frame(AUROC_scores) %>%
  rownames_to_column(var = "GO_term") %>%
  gather(Celltype, AUROC, -GO_term)
head(long_data) 
# GO_term   Celltype AUROC
# 1 GO:0000228 HTR3A.DRD2 0.915
# 2 GO:0000902 HTR3A.DRD2 0.996
# 3 GO:0000988 HTR3A.DRD2 0.964
# 4 GO:0003013 HTR3A.DRD2 0.991
# 5 GO:0003729 HTR3A.DRD2 0.848
# 6 GO:0003735 HTR3A.DRD2 0.657

ggplot(long_data, aes(x = Celltype, y = AUROC)) +
  geom_violin(aes(fill = Celltype), trim = FALSE) +  # Violin plot
  geom_jitter(width = 0.2, size = 2, alpha = 0.3) +  # Individual data points
  stat_summary(aes(group = interaction(Celltype)), fun = mean, geom = "crossbar", color = "black", width=.5, position = position_dodge(width = 0.8)) +
  labs(title = "AUROC scores of interneurons using GO terms", y = "AUROC", x = "Celltype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  coord_cartesian(ylim = c(0, 1))  # Set y-axis limits to 0 and 1

ggsave(here(plot_dir, "Metaneighbor_AUROC_interneuron_matched.pdf"), width=12, height=5)



# ===== Randomized list ===== 

sce.bla <- sce.backup

# Define the original "Inh_" cell types and their corresponding renaming
inh_celltypes <- c("Inh_5", "Inh_LAMP5_2", "Inh_4", "Inh_LAMP5_1", "Inh_ITC_2", "Inh_ITC_1", "Inh_2", "Inh_SST", "Inh_PVALB", "Inh_1", "Inh_VIP_1", "Inh_VIP_2", "Inh_3")
renaming <- c("PVALB ADAMTS5", "LAMP5 COL14A1", "HTR3A DRD2", "LAMP5 NDNF", "TSHZ1 SEMA3C", "TSHZ1 CALCRL", "VIP NDNF", "SST EPYC", "SST HGF", "VIP NDNF", "VIP ABI3BP", "VIP ABI3BP", "VIP ABI3BP")

# Randomize the renaming
set.seed(1234)  # Setting a seed for reproducibility
random_renaming <- sample(renaming)

# Apply the randomized renaming to the sce.bla object
for (i in seq_along(inh_celltypes)) {
  colData(sce.bla)$cell_type[colData(sce.bla)$cell_type == inh_celltypes[i]] <- random_renaming[i]
}

celltypes_of_interest = sce.bla[, sce.bla$cell_type %in% c('TSHZ1 SEMA3C', 'PVALB ADAMTS5', 'VIP ABI3BP',
                                                           'LAMP5 COL14A1', 'HTR3A DRD2', 'LAMP5 NDNF', 
                                                           'TSHZ1 CALCRL', 'VIP NDNF', 'SST EPYC', 
                                                           'SST HGF', 'VIP ABI3BP')]
celltype_matrix = model.matrix(~celltypes_of_interest$cell_type - 1)
colnames(celltype_matrix) = levels(as.factor(celltypes_of_interest$cell_type))

#data(GOhuman)
pdf(here(plot_dir,"Metaneighbor_AUROC_interneurons_shuffled.pdf"))
AUROC_scores_shuffled = MetaNeighbor(dat = celltypes_of_interest ,
                            experiment_labels = as.numeric(factor(celltypes_of_interest$study_id)),
                            celltype_labels = celltype_matrix,
                            genesets = GOhuman,
                            bplot = TRUE,
                            fast_version = TRUE)
dev.off()


# plotting using ggplot
long_data <- data.frame(AUROC_scores_shuffled) %>%
  rownames_to_column(var = "GO_term") %>%
  gather(Celltype, AUROC, -GO_term)
head(long_data) 
# GO_term   Celltype AUROC
# 1 GO:0000228 HTR3A.DRD2 0.915
# 2 GO:0000902 HTR3A.DRD2 0.996
# 3 GO:0000988 HTR3A.DRD2 0.964
# 4 GO:0003013 HTR3A.DRD2 0.991
# 5 GO:0003729 HTR3A.DRD2 0.848
# 6 GO:0003735 HTR3A.DRD2 0.657

ggplot(long_data, aes(x = Celltype, y = AUROC)) +
  geom_violin(aes(fill = Celltype), trim = FALSE) +  # Violin plot
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +  # Individual data points
  labs(title = "Violin plot of AUROC scores", y = "AUROC", x = "Celltype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  coord_cartesian(ylim = c(0, 1))  # Set y-axis limits to 0 and 1


# ============ Directly comparing real vs shuffled data ==============


# Convert AUROC_scores to long format and add a 'Type' column
long_real <- data.frame(AUROC_scores) %>%
  rownames_to_column(var = "GO_term") %>%
  gather(Celltype, AUROC, -GO_term) %>%
  mutate(Type = "Matched")

# Convert AUROC_scores_shuffled to long format and add a 'Type' column
long_shuffled <- data.frame(AUROC_scores_shuffled) %>%
  rownames_to_column(var = "GO_term") %>%
  gather(Celltype, AUROC, -GO_term) %>%
  mutate(Type = "Shuffled")

# Combine the two datasets
combined_data <- bind_rows(long_real, long_shuffled)

# Plot the combined data
ggplot(combined_data, aes(x = Celltype, y = AUROC)) +
  geom_violin(aes(fill = Type), trim = FALSE, position = position_dodge(width = 0.8), alpha=.3) +
  geom_jitter(aes(color = Type), position = position_dodge(width = 0.8), size = 2, alpha = 0.5) +
  stat_summary(aes(group = interaction(Type, Celltype)), fun = mean, geom = "crossbar", color = "black", width=.5, position = position_dodge(width = 0.8)) +
  labs(title = "Comparison of Matched vs Shuffled AUROC scores using GO terms", y = "AUROC", x = "Celltype") +
  scale_fill_manual(values = c("Matched" = "dodgerblue", "Shuffled" = "darkgrey")) +
  scale_color_manual(values = c("Matched" = "dodgerblue", "Shuffled" = "darkgrey")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 1))

ggsave(here(plot_dir, "Metaneighbor_AUROC_interneuron_matched_vs_shuffled.pdf"), width=12, height=5)


# =========== Functional characterization =================

library(GO.db)

# List of GO terms
go_terms <- c("GO:0030674", "GO:0030705", "GO:0043473")

# Retrieve GO term names
go_names <- select(GO.db, keys = as.character(names(GOhuman)), columns = "TERM", keytype = "GOID")

# Print the names
head(go_names)
# GOID                               TERM
# 1 GO:0000228                 nuclear chromosome
# 2 GO:0000902                 cell morphogenesis
# 3 GO:0000988                               <NA>
# 4 GO:0003013         circulatory system process
# 5 GO:0003729                       mRNA binding
# 6 GO:0003735 structural constituent of ribosome


gs_size = sapply(GOhuman, length)
aurocs_df = data.frame(go_term = go_names$TERM, AUROC_scores) 
aurocs_df$average = rowMeans(AUROC_scores)
aurocs_df$n_genes = gs_size[rownames(AUROC_scores)]

head(aurocs_df[order(aurocs_df$average, decreasing = TRUE),],10)
#                                                             go_term HTR3A.DRD2 LAMP5.COL14A1 LAMP5.NDNF
# GO:0000902                                       cell morphogenesis      0.996         0.998      0.995
# GO:0048646 anatomical structure formation involved in morphogenesis      0.993         0.999      0.995
# GO:0005578                                                     <NA>      0.988         0.998      0.993
# GO:0003013                               circulatory system process      0.991         0.995      0.985
# GO:0008283                            cell population proliferation      0.977         0.983      0.975
# GO:0040007                                                   growth      0.978         0.988      0.983
# GO:0008289                                            lipid binding      0.985         0.965      0.971
# GO:0030198                        extracellular matrix organization      0.968         0.996      0.988
# GO:0016757                             glycosyltransferase activity      0.963         0.984      0.978
# GO:0034330                               cell junction organization      0.967         0.974      0.949


small_sets = aurocs_df[aurocs_df$n_genes < 200,] 
head(small_sets[order(small_sets$average, decreasing = TRUE),],10)

#                                                   go_term HTR3A.DRD2 LAMP5.COL14A1 LAMP5.NDNF PVALB.ADAMTS5
# GO:0030674         protein-macromolecule adaptor activity      0.930         0.882      0.924         0.957
# GO:0030705 cytoskeleton-dependent intracellular transport      0.863         0.795      0.916         0.935
# GO:0043473                                   pigmentation      0.824         0.828      0.927         0.897
# GO:0007009                   plasma membrane organization      0.869         0.896      0.829         0.898
# GO:0042393                                histone binding      0.840         0.751      0.872         0.951
# GO:0008565                                           <NA>      0.859         0.840      0.833         0.861
# GO:0007034                             vacuolar transport      0.760         0.714      0.779         0.925
# GO:0032182                 ubiquitin-like protein binding      0.675         0.760      0.757         0.734
# GO:0016779                nucleotidyltransferase activity      0.718         0.756      0.882         0.783
# GO:0005811                                  lipid droplet      0.763         0.748      0.791         0.789
       