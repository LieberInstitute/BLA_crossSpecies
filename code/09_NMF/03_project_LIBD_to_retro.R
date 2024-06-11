library(RcppML)
library(here)
library(ggplot2)
library(projectR)
library(pheatmap)

plot_dir <- here("plots","09_NMF")
processed_dir <- here("processed-data", "09_NMF")

# ===== Load data =====
x <- readRDS(here(processed_dir, "RcppML_NMF_LIBD_Excit.rds"))

sce.human<- readRDS(here(processed_dir, "sce.excit.human.rds"))
sce.human

sce.retro<- readRDS(here(processed_dir, "sce.retro.excit.rds"))
sce.retro

counts(sce.retro) <- exp(assay(sce.retro))
logcounts(sce.retro) <- log(counts(sce.retro))


# check for duplicate colnames
any(duplicated(colnames(sce.retro)))

# make a hist of counts 0<x<1
hist(logcounts(sce.retro)[logcounts(sce.retro) > 0 & logcounts(sce.retro) < 1], breaks = 100)

# ===== Projecting patterns to retro seq =====
# == extract patterns ==
patterns <- t(x$h)
colnames(patterns) <- paste("NMF", 1:100, sep = "_")

loadings <- x$w
rownames(loadings) <- names(sce.human)

# == add to sce human ==

reducedDim(sce.human, "NMF_proj") <- patterns

for (i in 1:100){
    colData(sce.human)[[paste0("NMF_",i)]] <- reducedDims(sce.human)$NMF_proj[,i]
}


# == project  ==
rownames(sce.retro) <- rowData(sce.retro)$human_ortholog

sce.retro <- sce.retro[rownames(sce.retro) %in% rownames(loadings),]
sce.retro <- sce.retro[!duplicated(rownames(sce.retro)), ]
loadings <- loadings[rownames(loadings) %in% rownames(sce.retro),]
logcounts <- logcounts(sce.retro)
logcounts_mat <- Matrix(logcounts, sparse = TRUE)

proj <- project(logcounts_mat, loadings)
proj <- t(proj)
colnames(proj) <- paste("NMF", 1:100, sep = "_")

reducedDim(sce.retro, "NMF_proj") <- proj

# add each proj column to colData(spe)
for (i in 1:100){
    colData(sce.retro)[[paste0("NMF_",i)]] <- reducedDims(sce.retro)$NMF_proj[,i]
}

# drop duplicate colData
colData(sce.retro) <- colData(sce.retro)[ , !duplicated(colnames(colData(sce.retro)))]


# ======== Heatmaps ========

# == Human cell-types ==

# create dataframe 
data <- data.frame(colData(sce.human))

# aggregate NMF patterns across clusters. # grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("NMF", colnames(data))],
                      by=list(data$fine_type), 
                      FUN=mean)

# move Group.1 to row names, then drop
rownames(agg_data) <- agg_data$Group.1
agg_data <- agg_data[,-1]

p1 <- pheatmap(agg_data,
               color=colorRampPalette(c("blue","white","red"))(100),
               cluster_cols=T,
               cluster_rows=T,
               scale="column"
)

pdf(here(plot_dir, "Heatmap_NMF_human_patterns.pdf"), width=20, height=5)
p1
dev.off()


# == Exporting top gene per factor ==

# function for getting top n genes for each pattern
top_genes <- function(W, n=15){
    top_genes <- apply(W, 2, function(x) names(sort(x, decreasing=TRUE)[1:n]))
    return(top_genes)
}

#get top genes
top15 <- top_genes(loadings, 15)
head(top15)

# add colnames
colnames(top15) <- colnames(patterns)

# export to csv
write.csv(top15, here(processed_dir, "top15_genes_per_factor.csv"))


# ===== Retro seq heatmaps =====

# create dataframe 
data <- data.frame(colData(sce.retro)[,grep("NMF", colnames(data))])

# aggregate NMF patterns across clusters. # grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("NMF", colnames(data))],
                      by=list(data$Subclass), 
                      FUN=mean)

# move Group.1 to row names, then drop
rownames(agg_data) <- agg_data$Group.1
agg_data <- agg_data[,-1]


# drop any factos that sum 0
agg_data <- agg_data[, colSums(agg_data) != 0]

# normalize loadings
column_sums <- colSums(agg_data)

# Normalize each column by its sum
normalized_loadings <- sweep(agg_data, 2, column_sums, FUN="/")

p1 <- pheatmap(normalized_loadings,
               color=colorRampPalette(c("white","black"))(100),
               cluster_cols=T,
               cluster_rows=T
)

pdf(here(plot_dir, "Heatmap_NMF_retro_projection.pdf"), width=20, height=5)
p1
dev.off()



# ======= Manually selecting/annotating factors =======

# === Looking at all factors ===
nmf_df <- data.frame(colData(sce.human)[,grep("NMF", colnames(colData(sce.human)))])

# normalize loadings
col_sums <- colSums(nmf_df)

# Normalize each column by its sum
normalized_df <- sweep(nmf_df, 2, col_sums, FUN="/")
normalized_df$fine_type <- colData(sce.human)$fine_type
normalized_df <- normalized_df[!is.na(normalized_df) != 0,]

group <- normalized_df$fine_type
normalized_df$fine_type <- NULL

summarized <- summarizeAssayByGroup(t(normalized_df), 
                                    ids=DataFrame(group=group), 
                                    #subset.row=subset_patterns, 
                                    statistics=c("sum", "prop.detected"), threshold=0)

sum <- assay(summarized, "sum")
num <- assay(summarized, "prop.detected")
group.names <- summarized$group

evals_long <- data.frame(
    Feature=rep(colnames(normalized_df), ncol(num)),
    Group=rep(group.names, each=nrow(num)),
    NumDetected=as.numeric(num),
    Sum=as.numeric(sum)
)

p <- ggplot(evals_long) + 
    geom_point(aes_string(x="Group", y="Feature", size="NumDetected", col="Sum")) +
    scale_size(limits=c(0, max(evals_long$NumDetected))) + 
    scale_color_gradient(limits=c(0, max(evals_long$Sum)),
                         low="white", high="black") +
    theme(panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(size=0.5, colour = "grey80"),
          panel.grid.minor = element_line(size=0.25, colour = "grey80")) 

png(here(plot_dir, "DotPlot_NMF_human_patterns_vs_celltype.png"), width=15, height=30, units="in", res=300)
p
dev.off()

# === Selected patterns to use
# patterns to use:

# GPC5_ESR1: NMF_60
# GRM8_ADARB2: NMF_41
# GULP1_NTNG1: NMF_9
# CHRM3_TRPS1: NMF_30
# COL25A1_PEX5L: NMF_40
# GULP1_ZBTB20: NMF_49


# subset to these factors in agg_data
subset_patterns <- c("NMF_60", "NMF_41", "NMF_9", "NMF_30", "NMF_40", "NMF_49")



nmf_df <- data.frame(colData(sce.human)[,grep("NMF", colnames(colData(sce.human)))])



# normalize loadings
col_sums <- colSums(nmf_df)

# Normalize each column by its sum
normalized_df <- sweep(nmf_df, 2, col_sums, FUN="/")
normalized_df$fine_type <- colData(sce.human)$fine_type
normalized_df <- normalized_df[!is.na(normalized_df) != 0,]

# drop ZNF804B_TTN, TENM3_SATB2, and SAMD5_NTNG1 celltypes
normalized_df <- normalized_df[!normalized_df$fine_type %in% c("ZNF804B_TTN", "TENM3_SATB2", "SAMD5_NTNG1"),]

group <- normalized_df$fine_type
normalized_df$fine_type <- NULL



summarized <- summarizeAssayByGroup(t(normalized_df), 
                                    ids=DataFrame(group=group), 
                                    subset.row=subset_patterns, 
                                    statistics=c("sum", "prop.detected"), threshold=0)

sum <- assay(summarized, "sum")
num <- assay(summarized, "prop.detected")
group.names <- summarized$group

evals_long <- data.frame(
    Feature=rep(subset_patterns, ncol(num)),
    Group=rep(group.names, each=nrow(num)),
    NumDetected=as.numeric(num),
    Sum=as.numeric(sum)
)

#re-order factors for better plottin
evals_long$Feature <- factor(evals_long$Feature, levels = c("NMF_49", "NMF_9", "NMF_41", "NMF_60", "NMF_40", "NMF_30"))

p <- ggplot(evals_long) + 
    geom_point(aes_string(x="Group", y="Feature", size="NumDetected", col="Sum")) +
    scale_size(name="proportion\nnuclei withn\nnon-zero\nweight",
               limits=c(0, max(evals_long$NumDetected))) + 
    scale_color_gradient(name="aggregate\npredicted\nnuclei-level\nweights",
                         limits=c(0, max(evals_long$Sum)),
                         low="white", high="black") +
    theme(panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(size=0.5, colour = "grey80"),
          panel.grid.minor = element_line(size=0.25, colour = "grey80")) +
    coord_flip() + 
    xlab("Human Cell-types") + ylab("NMF pattern")


png(here(plot_dir, "DotPlot_NMF_human_subset_patterns.png"), width=8, height=4.5, units="in", res=300)
p
dev.off()


# ========= Dot plots of NMF loadings vs retro seq subclass ==========
nmf_df <- data.frame(colData(sce.retro)[,grep("NMF", colnames(colData(sce.retro)))])

# normalize loadings
col_sums <- colSums(nmf_df)

# Normalize each column by its sum
normalized_df <- sweep(nmf_df, 2, col_sums, FUN="/")
normalized_df$Subclass <- colData(sce.retro)$Subclass
normalized_df <- normalized_df[!is.na(normalized_df) != 0,]

group <- normalized_df$Subclass
normalized_df$Subclass <- NULL

summarized <- summarizeAssayByGroup(t(normalized_df), 
                                    ids=DataFrame(group=group), 
                                    subset.row=subset_patterns, 
                                    statistics=c("sum", "prop.detected"), threshold=0)

sum <- assay(summarized, "sum")
num <- assay(summarized, "prop.detected")
group.names <- summarized$group

evals_long <- data.frame(
    Feature=rep(subset_patterns, ncol(num)),
    Group=rep(group.names, each=nrow(num)),
    NumDetected=as.numeric(num),
    Sum=as.numeric(sum)
)

#evals_long$Group <- factor(evals_long$Group, levels = c(" L2/3 IT PIR-ENTl Glut", "MEA Slc17a7 Glut","NLOT Rho Glut", 
#                                                        "SI-MA-LPO-LHA Skor1 Glut", "MEA-BST Otp Zic2 Glut",
#                                                        "COAa-PAA-MEA Barhl2 Glut","MEA-COA-BMA Codd42 Glut","MEA Otp Foxp2 Glut",
#                                                        "PVH-SO-PVa Otp Glut"
#                                                        ))

p<-ggplot(evals_long) + 
    geom_point(aes_string(x="Group", y="Feature", size="NumDetected", col="Sum")) +
    scale_size(name="proportion\nnuclei withn\nnon-zero\nweight",
               limits=c(0, max(evals_long$NumDetected))) + 
    scale_color_gradient(name="aggregate\npredicted\nnuclei-level\nweights",
                         limits=c(0, max(evals_long$Sum)),
                         low="white", high="black") +
    theme(panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(size=0.5, colour = "grey80"),
          panel.grid.minor = element_line(size=0.25, colour = "grey80")) +
    coord_flip() + 
    xlab("Allen Institute Subclass") + ylab("NMF pattern")

png(here(plot_dir, "DotPlot_NMF_retro_patterns_vs_celltype.png"), width=8, height=5.5, units="in", res=300)
p
dev.off()


# ===== Target region =====
nmf_df <- data.frame(colData(sce.retro)[,grep("NMF", colnames(colData(sce.retro)))])

# normalize loadings
col_sums <- colSums(nmf_df)

# Normalize each column by its sum
normalized_df <- sweep(nmf_df, 2, col_sums, FUN="/")
normalized_df$Subclass <- colData(sce.retro)$Target
normalized_df <- normalized_df[!is.na(normalized_df) != 0,]

group <- normalized_df$Subclass
normalized_df$Subclass <- NULL

summarized <- summarizeAssayByGroup(t(normalized_df), 
                                    ids=DataFrame(group=group), 
                                    #subset.row=subset_patterns, 
                                    statistics=c("sum", "prop.detected"), threshold=0)

sum <- assay(summarized, "sum")
num <- assay(summarized, "prop.detected")
group.names <- summarized$group

evals_long$Group <- factor(evals_long$Group, levels = c())

evals_long <- data.frame(
    Feature=rep(colnames(normalized_df), ncol(num)),
    Group=rep(group.names, each=nrow(num)),
    NumDetected=as.numeric(num),
    Sum=as.numeric(sum)
)

#re-order factors for better plottin
evals_long$Feature <- factor(evals_long$Feature, levels = c("NMF_49", "NMF_9", "NMF_41", "NMF_60", "NMF_40", "NMF_30"))

p <- ggplot(evals_long) + 
    geom_point(aes_string(x="Group", y="Feature", size="NumDetected", col="Sum")) +
    scale_size(name="proportion\nnuclei withn\nnon-zero\nweight",
               limits=c(0, max(evals_long$NumDetected))) + 
    scale_color_gradient(name="aggregate\npredicted\nnuclei-level\nweights",
                         limits=c(0, max(evals_long$Sum)),
                         low="white", high="black") +
    theme(panel.background = element_rect(fill="white"),
          panel.grid.major = element_line(size=0.5, colour = "grey80"),
          panel.grid.minor = element_line(size=0.25, colour = "grey80")) +
    #coord_flip() +
    xlab("Target Region") + ylab("NMF pattern")

png(here(plot_dir, "DotPlot_NMF_retro_patterns_vs_target_all.png"), width=5, height=30, units="in", res=300)
p
dev.off()
