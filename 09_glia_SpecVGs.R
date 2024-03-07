devtools::install("/Users/mtotty2/Documents/R/SpecVGs")
library(SpecVGs)

library(STexampleData)
library(SpatialExperiment)
library(here)
library(ggplot2)
library(scuttle)
library(scater)
library(scran)

# NOTE:: This failed due to exhausted memory. lm requires a dense matrix, but it looks like
# MatrixModels::lm.fit.sparse() can run lm on sparse matrices. Should look into this
# as next steps.
#
# Ref: https://stackoverflow.com/questions/3169371/large-scale-regression-in-r-with-a-sparse-feature-matrix


# load data



# ====== Formatting for SpecVGs ======


# Example assuming 'Species' and 'Organ' are column names in your object
gene_names <- rownames(assay(sce))  # Get gene names
results <- list()  # To store results for each gene

# make data frame with genes as rows, species as columns, with logcounts
# as values
mat <- as.matrix(logcounts(sce))
df <- as.data.frame(t(mat))

# get species and celltype data
species_col <- colData(sce)$species
celltype <- colData(sce)$broad_celltype
ident <- colData(sce)$ident

# add to dataframe
df$species <- species_col
df$celltype <- celltype
df$ident <- ident

# Replace "-" with "_" in gene names
gene_names <- gsub("-", "_", gene_names)
colnames(df) <- gsub("-", "_", colnames(df))


# ===== get gene variances =====
startTime <- Sys.time()
var_ident <- SpecVGs::calculate_gene_variances(gene_names=gene_names,
                                               df=df,
                                               covariates = c("ident", "species"),
                                               verbose=TRUE
)
endTime <- Sys.time()

# prints recorded time 
print(endTime - startTime)

# calc total
var_ident$total <- var_ident$species + var_ident$ident + var_ident$Resid_Variance

# do this with var_ident
var_ident$species_norm <- var_ident$species / var_ident$total
var_ident$ident_norm <- var_ident$ident / var_ident$total

head(var_ident)
# Gene Resid_Variance     ident     species     total species_norm   ident_nom
# 1  ANKRD65       461.5926 10.042796  3.96300490  475.5984 0.0083326704 0.021116127
# 2 AURKAIP1       686.3674 11.757966 10.95959978  709.0849 0.0154559760 0.016581887
# 3    MXRA8       251.4145  6.017750  0.04264901  257.4748 0.0001656434 0.023372184
# 4     DVL1      2456.2776 31.383257  6.28085157 2493.9417 0.0025184437 0.012583798
# 5   TAS1R3       524.4282  3.369079  5.24777505  533.0451 0.0098448994 0.006320440
# 6     CPTP       291.5154  1.902968  3.16581772  296.5842 0.0106742618 0.006416283

# ====== plotting ======

plot(var_ident$species_norm, var_ident$ident_norm, 
     xlab = "Species Variance", ylab = "Celltype Variance", 
     main = "Scatter Plot of Species vs Celltype Variance",
     pch = 19, xlim = c(0,1), ylim = c(0,1), cex=.5)


abline(coef = c(0.5, -1), col = "grey", lwd = 2)
abline(coef = c(0, 2), col = "grey", lwd = 2)
abline(coef = c(0, 0.5), col = "grey", lwd = 2)



# ========= Selecting SVGs vs TVGs ===========
# SVGs (species variable genes) should be defined as celltype + species variance > 0.5
# AND species > 2*celltype

# TVGs (total variable genes) should be defined as celltype + species variance > 0.5
# AND species < 2*celltype

# Define SVGs
svg_genes <- var_ident[which( (var_ident$species_norm + var_ident$ident_norm) > 0.5 & var_ident$species_norm > 2*var_ident$ident_norm), ]
dim(svg_genes)[1]
# [1] 20

# Define TVGs
tvg_genes <- var_ident[which( (var_ident$species_norm + var_ident$ident_norm) > 0.5 & var_ident$ident_norm > 2*var_ident$species_norm), ]
dim(tvg_genes)[1]
# [1] 303

points(tvg_genes$species_norm, tvg_genes$ident_norm, pch = 16, col = "orange", cex=.6)
points(svg_genes$species_norm, svg_genes$ident_norm, pch = 16, col = "darkgreen", cex=.6)
