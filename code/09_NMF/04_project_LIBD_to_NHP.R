library(RcppML)
library(here)
library(ggplot2)
library(projectR)
library(pheatmap)
library(Matrix)
library(scater)
library(SummarizedExperiment)

plot_dir <- here("plots","09_NMF")
processed_dir <- here("processed-data", "09_NMF")

# ===== Load data =====
x <- readRDS(here(processed_dir, "RcppML_NMF_LIBD_Excit.rds"))

sce.human<- readRDS(here(processed_dir, "sce.excit.human.rds"))
sce.human

sce.excit <- readRDS(here("processed-data", "07_annotation", "sce.excit.integrated.annotated.rds"))
sce.baboon <- sce.excit[, sce.excit$species == "baboon"]
sce.baboon

sce.mac <- sce.excit[, sce.excit$species == "macaque"]
sce.mac


# ===== Projecting patterns to mac seq =====
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


# == project to macaque  ==
#rownames(sce.mac) <- rowData(sce.mac)$human_ortholog

sce.mac <- sce.mac[rownames(sce.mac) %in% rownames(loadings),]
sce.mac <- sce.mac[!duplicated(rownames(sce.mac)), ]
loadings <- loadings[rownames(loadings) %in% rownames(sce.mac),]
logcounts <- logcounts(sce.mac)
logcounts_mat <- Matrix(logcounts, sparse = TRUE)

proj <- project(logcounts_mat, loadings)
proj <- t(proj)
colnames(proj) <- paste("NMF", 1:100, sep = "_")

reducedDim(sce.mac, "NMF_proj") <- proj

# add each proj column to colData(spe)
for (i in 1:100){
    colData(sce.mac)[[paste0("NMF_",i)]] <- reducedDims(sce.mac)$NMF_proj[,i]
}

# drop duplicate colData
colData(sce.mac) <- colData(sce.mac)[ , !duplicated(colnames(colData(sce.mac)))]



# == project to baboon  ==
sce.baboon <- sce.baboon[rownames(sce.baboon) %in% rownames(loadings),]
sce.baboon <- sce.baboon[!duplicated(rownames(sce.baboon)), ]
loadings <- loadings[rownames(loadings) %in% rownames(sce.baboon),]
logcounts <- logcounts(sce.baboon)
logcounts_mat <- Matrix(logcounts, sparse = TRUE)

proj <- project(logcounts_mat, loadings)
proj <- t(proj)
colnames(proj) <- paste("NMF", 1:100, sep = "_")

reducedDim(sce.baboon, "NMF_proj") <- proj

# add each proj column to colData(spe)
for (i in 1:100){
    colData(sce.baboon)[[paste0("NMF_",i)]] <- reducedDims(sce.baboon)$NMF_proj[,i]
}

# drop duplicate colData
colData(sce.baboon) <- colData(sce.baboon)[ , !duplicated(colnames(colData(sce.baboon)))]


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


# ===== mac seq heatmaps =====

# create dataframe 
data <- data.frame(colData(sce.mac))

# aggregate NMF patterns across clusters. # grep "NMF" to get all NMF patterns
agg_data <- aggregate(data[,grep("NMF", colnames(data))],
                      by=list(data$fine_type), 
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

pdf(here(plot_dir, "Heatmap_NMF_mac_projection.pdf"), width=20, height=5)
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


# ========= Dot plots of NMF loadings vs mac seq subclass ==========
nmf_df <- data.frame(colData(sce.mac)[,grep("NMF", colnames(colData(sce.mac)))])

# normalize loadings
col_sums <- colSums(nmf_df)

# Normalize each column by its sum
normalized_df <- sweep(nmf_df, 2, col_sums, FUN="/")
normalized_df$fine_type <- colData(sce.mac)$fine_type
normalized_df <- normalized_df[!is.na(normalized_df) != 0,]

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
evals_long$Feature <- factor(evals_long$Feature, 
                             levels = c("NMF_49", "NMF_9", "NMF_41", "NMF_60", "NMF_40", "NMF_30")
                             )


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

png(here(plot_dir, "DotPlot_NMF_mac_patterns_vs_celltype.png"), width=8, height=5.5, units="in", res=300)
p
dev.off()


# ========= Dot plots of NMF loadings vs baboon seq subclass ==========
nmf_df <- data.frame(colData(sce.baboon)[,grep("NMF", colnames(colData(sce.baboon)))])

# normalize loadings
col_sums <- colSums(nmf_df)

# Normalize each column by its sum
normalized_df <- sweep(nmf_df, 2, col_sums, FUN="/")
normalized_df$fine_type <- colData(sce.baboon)$fine_type
normalized_df <- normalized_df[!is.na(normalized_df) != 0,]

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
evals_long$Feature <- factor(evals_long$Feature, 
                             levels = c("NMF_49", "NMF_9", "NMF_41", "NMF_60", "NMF_40", "NMF_30")
)


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

png(here(plot_dir, "DotPlot_NMF_baboon_patterns_vs_celltype.png"), width=8, height=5.5, units="in", res=300)
p
dev.off()