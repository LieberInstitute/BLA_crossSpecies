library(dplyr)
library(here)
library(patchwork)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(EnhancedVolcano)
library(tidyr)
library(scater)
library(Seurat)


# directories
processed_dir = here("processed-data", "08_species_comparisons")
plot_dir = here("plots", "08_species_comparisons", "03_DE_analysis")

biccn <- readRDS(here(processed_dir, "BICCN_Amy_Excit_Supercluster.rds"))
biccn

# convert to SCE
sce <- as.SingleCellExperiment(biccn)
sce
# class: SingleCellExperiment 
# dim: 59357 109452 
# metadata(0):
#     assays(1): logcounts
# rownames(59357): ENSG00000256637 ENSG00000118523 ... ENSG00000164180 ENSG00000273449
# rowData names(0):
#     colnames(109452): 10X349_4:CCCATTGTCCGCTTAC 10X349_3:TTCTCTCTCTCCCATG ... 10X203_3:ATATCCTTCAGTAGGG
# 10X270_8:TGGGAGACAGGATCTT
# colData names(33): roi organism_ontology_term_id ... observation_joinid ident
# reducedDimNames(2): UMAP TSNE
# mainExpName: RNA
# altExpNames(0):

colnames(colData(sce))
# [1] "roi"                                      "organism_ontology_term_id"               
# [3] "disease_ontology_term_id"                 "self_reported_ethnicity_ontology_term_id"
# [5] "assay_ontology_term_id"                   "sex_ontology_term_id"                    
# [7] "development_stage_ontology_term_id"       "donor_id"                                
# [9] "suspension_type"                          "dissection"                              
# [11] "fraction_mitochondrial"                   "fraction_unspliced"                      
# [13] "cell_cycle_score"                         "total_genes"                             
# [15] "total_UMIs"                               "sample_id"                               
# [17] "supercluster_term"                        "cluster_id"                              
# [19] "subcluster_id"                            "cell_type_ontology_term_id"              
# [21] "tissue_ontology_term_id"                  "is_primary_data"                         
# [23] "tissue_type"                              "cell_type"                               
# [25] "assay"                                    "disease"                                 
# [27] "organism"                                 "sex"                                     
# [29] "tissue"                                   "self_reported_ethnicity"                 
# [31] "development_stage"                        "observation_joinid"                      
# [33] "ident"       


sce$dissection <- droplevels(sce$dissection)
unique(sce$dissection)


# ==== plot UMAP

plotReducedDim(sce, dimred="UMAP", colour_by="dissection") +
    theme(legend.position="none")


# get only dissection with "AMY" in them
sce.amy <- sce[, grepl("AMY", sce$dissection)]
sce.amy
# class: SingleCellExperiment 
# dim: 59357 61749 
# metadata(0):
#     assays(1): logcounts
# rownames(59357): ENSG00000256637 ENSG00000118523 ... ENSG00000164180 ENSG00000273449
# rowData names(0):
#     colnames(61749): 10X191_7:ACTGATGAGAAGAACG 10X350_7:AAAGGGCCAACCCGCA ... 10X203_5:TTATTGCGTACCCGAC
# 10X203_5:TGAGCGCAGGGAGGGT
# colData names(33): roi organism_ontology_term_id ... observation_joinid ident
# reducedDimNames(2): UMAP TSNE
# mainExpName: RNA
# altExpNames(0):

sce.amy$dissection <- droplevels(sce.amy$dissection)
unique(sce.amy$dissection)

plotReducedDim(sce.amy, dimred="UMAP", colour_by="dissection") +
    theme(legend.position="none")


dissection_labels <- unique(sce.amy$dissection)
# [1] Amygdaloid complex (AMY) - Corticomedial nuclear group (CMN) - anterior cortical nucleus - CoA                
# [2] Amygdaloid complex (AMY) - Central nuclear group - CEN                                                        
# [3] Amygdaloid complex (AMY) - corticomedial nuclear group - CMN                                                  
# [4] Amygdaloid complex (AMY) - basolateral nuclear group (BLN) - basolateral nucleus (basal nucleus) - BL         
# [5] Amygdaloid complex (AMY) - basolateral nuclear group (BLN) - basomedial nucleus (accessory basal nucleus) - BM
# [6] Amygdaloid complex (AMY) - Basolateral nuclear group (BLN) - lateral nucleus - La                             
# 6 Levels: Amygdaloid complex (AMY) - Basolateral nuclear group (BLN) - lateral nucleus - La ...

# name the dissecton levels accordingly
# [1] <- CoA
# [2] <- CeA
# [3] <- CoM
# [4] <- BL
# [5] <- BM
# [6] <- LA
sce.amy$subregion <- c()

sce.amy$subregion[sce.amy$dissection == dissection_labels[1]] <- "CoA"
sce.amy$subregion[sce.amy$dissection == dissection_labels[2]] <- "CeA"
sce.amy$subregion[sce.amy$dissection == dissection_labels[3]] <- "CoM"
sce.amy$subregion[sce.amy$dissection == dissection_labels[4]] <- "BL"
sce.amy$subregion[sce.amy$dissection == dissection_labels[5]] <- "BM"
sce.amy$subregion[sce.amy$dissection == dissection_labels[6]] <- "LA"

plotReducedDim(sce.amy, dimred="UMAP", colour_by="subregion")


# ===== Convert from Ensembl ID to gene name =====
library("biomaRt")


# Step 1: Connect to the Ensembl BioMart database
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Step 2: Retrieve gene names and biotypes for your Ensembl IDs
genes <- rownames(sce.amy)
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
                   filters = 'ensembl_gene_id', values = genes, mart = mart)

# Step 3: Filter for protein-coding genes
protein_coding_genes <- gene_info[gene_info$gene_biotype == 'protein_coding',]

# Step 4: Filter the SCE object to keep only protein-coding genes
sce.amy <- sce.amy[protein_coding_genes$ensembl_gene_id, ]

# Step 5: Update row names with gene symbols
rowData(sce.amy)$Ensembl_ID <- rownames(sce.amy)
rowData(sce.amy)$Gene_name <- protein_coding_genes$external_gene_name

rownames(sce.amy) <- rowData(sce.amy)$Gene_name
sce.amy

saveRDS(sce.amy, here(processed_dir, "BICCN_processed.rds"))


# === rerun PCA and UMAP to clean up this data ===
sce.amy <- readRDS(here(processed_dir, "BICCN_processed.rds"))

# remove blank row names
sce.amy <- sce.amy[-which(rownames(sce.amy) == "",),]

# rename logcounts to counts
counts(sce.amy) <- logcounts(sce.amy)

# lognormalize
sce.amy <- logNormCounts(sce.amy)

# get HVGs
model <- modelGeneVar(sce.amy)
HVGs <- getTopHVGs(model, n=1000)

# run PCA
sce.amy <- scran::fixedPCA(sce.amy, subset.row=HVGs)


# batch correct 
# --- NOTES---
# need to look into multi-factor batch correction for both
# dissection and donor_id
test <- batchelor::fastMNN(sce.amy, batch=sce.amy$dissection)
reducedDims(sce.amy)$corrected <- reducedDims(test)$corrected

# run UMAP
sce.amy <- runUMAP(sce.amy, dimred="corrected", min_dist=0.3)

# plot new UMAP

p1 <- plotReducedDim(sce.amy, dimred="UMAP", colour_by="subregion", point_size=1) 
p2 <- plotReducedDim(sce.amy, dimred="UMAP", colour_by="GULP1", point_size=1) +
    scale_colour_gradient(low="grey", high="red") 
p3 <- plotReducedDim(sce.amy, dimred="UMAP", colour_by="COL25A1", point_size=1)  +
    scale_colour_gradient(low="grey", high="red")

(p1 + p2 + p3)


plotExpression(sce.amy, features=c("GULP1", "COL25A1","GAD1"), x="subregion", colour_by="subregion",
               ncol=1)

