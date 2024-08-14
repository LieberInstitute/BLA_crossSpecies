'''
Notes:

I currently have run this using the combined SCE object with all species. 
Need to subset to only humans.

ALSO - I should see if I can add back in all the genes dropped when combining to
1-1 orthologs. For example, one of the top DEGs CORT is not in NHP genes.

ALSO - After subsetting the human only Im getting some now weird results, 
which makes me think I need to re-normalize the human data.

'''



library("spatialLIBD")
library("patchwork")
library("scran")
library("fgsea")
library("here")
library("LRcell")
library(scater)

# Save directories
plot_dir = here("plots", "09_psychiatric_GSEA","LRcell")
processed_dir <- here("processed-data",  "09_psychiatric_GSEA")

# ====== Load SCE ======
sce <- readRDS(here("processed-data", "09_psychiatric_GSEA", "sce.human_all_genes.rds"))

# filter cell types with < 70 cells
celltype_counts <- table(colData(sce)$fine_celltype)
celltypes_to_drop <- names(celltype_counts[celltype_counts < 70])
sce <- sce[, !(colData(sce)$fine_celltype %in% celltypes_to_drop)]
sce$fine_celltype <- factor(sce$fine_celltype)


# ====== get Jaffe DEGs =====
jaffe_degs <- read.csv(here(processed_dir,"Jaffe_DEGs.csv"))



# =========================
#  PTSD
# =========================

# ===== Expression Vector =====
# Filter the significant genes (p-value less than 0.05)
PTSD_BLA_genes <- jaffe_degs[jaffe_degs$BLA_PValue_PTSD < 0.05,]

# Prepare the named vector for LRcell
PTSD_degs_vector <- setNames(PTSD_BLA_genes$BLA_PValue_PTSD, PTSD_BLA_genes$Symbol)
length(PTSD_degs_vector)
#[1] 1876

# Subset to only genes we have in the sce object
existing_genes <- rownames(sce)[rownames(sce) %in% names(PTSD_degs_vector)]

# Extract logcounts as matrix
expression_matrix <- as.matrix(logcounts(sce)[existing_genes,])



# ===== Annotation vector =====
# Extract the cell type annotations for each fine cell type
cell_type_annotations <- sce$fine_celltype

# Add cell IDs to the annotation vector
names(cell_type_annotations) <- colnames(sce)


# ==== LRcell enrichment scores =====
# generating the enrichment score 
enriched_res <- LRcell_gene_enriched_scores(expr = expression_matrix,
                                            annot = cell_type_annotations, parallel = FALSE)
enriched_res

# Convert to a dataframe, replaced NA with 0
df <- as.data.frame(t(enriched_res))
df[is.na(df)] <- 0

# Now calculate average of each row (which represent cell type)
df$avg <- rowMeans(df)


# ====== Running LRcell =======
# get marker genes for LRcell in logistic regression
BLA_marker_genes <- get_markergenes(enriched_res, method="LR", topn=25)

# view some of the marke genes for each cluster
head(lapply(BLA_marker_genes, head))

# run LRcell
PTSD_res <- LRcellCore(gene.p = PTSD_degs_vector,
                  marker.g = BLA_marker_genes,
                  method = "LR")

# clean up output for better visualization
PTSD_res$cell_type <- unlist(lapply(strsplit(PTSD_res$ID, '\\.'), '[', 2))
PTSD_res$cell_type <- unique(sce$fine_celltype) 





# =========================
#  MDD
# =========================
# Filter the significant genes (p-value less than 0.05)
MDD_BLA_genes <- jaffe_degs[jaffe_degs$BLA_PValue_MDD < 0.05,]

# Prepare the named vector for GSEA 
MDD_degs_vector <- setNames(MDD_BLA_genes$BLA_PValue_MDD, MDD_BLA_genes$Symbol)
length(MDD_degs_vector)
# [1] 1647

# order genes by p-value, descending
MDD_degs_vector <- MDD_degs_vector[order(MDD_degs_vector, decreasing = FALSE)]

# Subset to only genes we have in the sce object
existing_genes <- rownames(sce)[rownames(sce) %in% names(MDD_degs_vector)]

# Extract logcounts as matrix
expression_matrix <- as.matrix(logcounts(sce)[existing_genes,])


# ==== LRcell enrichment scores =====
# generating the enrichment score 
enriched_res <- LRcell_gene_enriched_scores(expr = expression_matrix,
                                            annot = cell_type_annotations, parallel = TRUE)
enriched_res

# Convert to a dataframe, replaced NA with 0
df <- as.data.frame(t(enriched_res))
df[is.na(df)] <- 0

# Now calculate average of each row (which represent cell type)
df$avg <- rowMeans(df)

# ====== Running LRcell =======
# get marker genes for LRcell in logistic regression
BLA_marker_genes <- get_markergenes(enriched_res, method="LR", topn=25)

# view some of the marke genes for each cluster
head(lapply(BLA_marker_genes, head))

# run LRcell
MDD_res <- LRcellCore(gene.p = MDD_degs_vector,
                  marker.g = BLA_marker_genes,
                  method = "LR")

# clean up output for better visualization
MDD_res$cell_type <- unlist(lapply(strsplit(MDD_res$ID, '\\.'), '[', 2))
MDD_res$cell_type <- unique(sce$fine_celltype) # add broad cell types here for x-axis




# ====== Plotting ======
# Code taken from LRcell::plot_manhattan_enrich function

label.topn <- 5
sig.cutoff <- .05

# PTSD
res <- PTSD_res
lr_tmp <- res %>%
        dplyr::group_by(.data$cell_type) %>%
        dplyr::arrange(.data$cell_type)
    lr_tmp$pos <- as.numeric(rownames(lr_tmp))

    # axis center
    axisdf <- lr_tmp %>%
        dplyr::group_by(.data$cell_type) %>%
        dplyr::summarize(center = (max(.data$pos) + min(.data$pos)) / 2)

    thres_topn <- ifelse(nrow(lr_tmp) > label.topn, label.topn, nrow(lr_tmp))
    repel_thres <- sort(lr_tmp$FDR)[thres_topn]

 p1 <- ggplot2::ggplot(lr_tmp, aes(x = -log10(.data$FDR), y = .data$pos)) +
    # add points
    ggplot2::geom_point(aes(color = as.factor(.data$cell_type)), alpha = 0.8, size = 2) +
    ggplot2::geom_vline(xintercept = -log10(sig.cutoff), linetype = "dashed", color = "red") +
    ggrepel::geom_text_repel(aes(label = ifelse(.data$FDR <= repel_thres & .data$FDR <= sig.cutoff, as.character(.data$ID), "")),
                    force = 5) +
    ggplot2::scale_y_continuous(label = axisdf$cell_type, breaks = axisdf$center) +
    # add labs
    ggplot2::labs(y = "Human Cell Types",
    title="PTSD DEGs from Jaffe et al. (2021)") +
    ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")


# MDD
res <- MDD_res
lr_tmp <- res %>%
        dplyr::group_by(.data$cell_type) %>%
        dplyr::arrange(.data$cell_type)
    lr_tmp$pos <- as.numeric(rownames(lr_tmp))

    # axis center
    axisdf <- lr_tmp %>%
        dplyr::group_by(.data$cell_type) %>%
        dplyr::summarize(center = (max(.data$pos) + min(.data$pos)) / 2)

    thres_topn <- ifelse(nrow(lr_tmp) > label.topn, label.topn, nrow(lr_tmp))
    repel_thres <- sort(lr_tmp$FDR)[thres_topn]

 p2 <- ggplot2::ggplot(lr_tmp, aes(x = -log10(.data$FDR), y = .data$pos)) +
    # add points
    ggplot2::geom_point(aes(color = as.factor(.data$cell_type)), alpha = 0.8, size = 2) +
    ggplot2::geom_vline(xintercept = -log10(sig.cutoff), linetype = "dashed", color = "red") +
    ggrepel::geom_text_repel(aes(label = ifelse(.data$FDR <= repel_thres & .data$FDR <= sig.cutoff, as.character(.data$ID), "")),
                    force = 5) +
    ggplot2::scale_y_continuous(label = axisdf$cell_type, breaks = axisdf$center) +
    # add labs
    ggplot2::labs(y = "Human Cell Types",
    title="MDD DEGs from Jaffe et al. (2021)") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
     axis.text.y=element_blank(),
     axis.ticks.y=element_blank(),
     axis.title.y=element_blank()
    )



png(here(plot_dir, "LRcell_manhattan.png"), width = 10, height = 10, units = "in", res = 300)
p1 + p2
dev.off()