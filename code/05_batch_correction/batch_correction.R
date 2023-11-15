library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)
library(batchelor)




## save directories
plot_dir = here("plots", "05_batch_correction")
processed_dir = here("processed-data","05_batch_correction")

# load sce
# save combined, uncorrected sce
load(here("processed-data","04_norm_and_dim_reduction", "sce_combined_uncorrected.rda"))
combined
# class: SingleCellExperiment 
# dim: 14391 187974 
# metadata(0):
#     assays(4): merged counts logcounts binomial_deviance_residuals
# rownames(14391): SAMD11 NOC2L ... EIF1AY RPS4Y2
# rowData names(4): gene_name gene_type Symbol.uniq gene_biotype
# colnames(187974): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# colData names(15): batch Sample ... doubletScore sizeFactor
# reducedDimNames(4): PCA TSNE GLM-PCA_approx UMAP
# mainExpName: NULL
# altExpNames(0):

unique(combined$batch)
# [1] "Human"   "Macaque" "Baboon" 

combined$Species <- combined$batch
unique(combined$Species)
# [1] "Human"   "Macaque" "Baboon" 


# ================= Batch correction with MNN =================

# Run batch correction with MNN across Samples
mnn_samples <- batchelor::reducedMNN(reducedDim(combined, "GLM-PCA_approx"),
                                    batch=as.factor(combined$Sample))

reducedDim(combined,"mnn_samples") <- mnn_samples$corrected

# Run batch correction with MNN across Species 
mnn_species <- batchelor::reducedMNN(reducedDim(combined, "GLM-PCA_approx"),
                                     batch=as.factor(combined$Species))

reducedDim(combined,"mnn_species") <- mnn_species$corrected

combined
# class: SingleCellExperiment 
# dim: 14391 187974 
# metadata(0):
#     assays(4): merged counts logcounts binomial_deviance_residuals
# rownames(14391): SAMD11 NOC2L ... EIF1AY RPS4Y2
# rowData names(4): gene_name gene_type Symbol.uniq gene_biotype
# colnames(187974): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# colData names(16): batch Sample ... sizeFactor Species
# reducedDimNames(6): PCA TSNE ... mnn_samples mnn_species
# mainExpName: NULL
# altExpNames(0):



# ============= TSNE and UMAP on MNN corrected data ===============

# calculate UMAP and TSNE using GLM-PCA_approx
set.seed(1234)
# combined <- runUMAP(combined,
#                     dimred = "mnn_species",
#                     name = "UMAP_mnn_species",
#                     #BPPARAM = BiocParallel::MulticoreParam(workers=20)
#                     )

combined <- runTSNE(combined,
                    dimred = "mnn_species",
                    name = "TSNE_mnn_species"
                    #BPPARAM = BiocParallel::MulticoreParam(workers=20)
                    )

pdf(here(plot_dir, "TSNE_mnn_corrected_species.pdf"))
plotReducedDim(combined, dim="TSNE_mnn_species", colour_by = "Species", point_alpha = 0.2)
dev.off()



combined <- runTSNE(combined,
                    dimred = "mnn_samples",
                    name = "TSNE_mnn_samples"
                    #BPPARAM = BiocParallel::MulticoreParam(workers=20)
)

pdf(here(plot_dir, "TSNE_mnn_corrected_samples.pdf"))
plotReducedDim(combined, dim="TSNE_mnn_samples", colour_by = "Species", point_alpha = 0.2)
dev.off()





# Run batch correction with MNN across Samples
mnn_samples <- batchelor::reducedMNN(reducedDim(combined, "PCA"),
                                     batch=as.factor(combined$Sample))

reducedDim(combined,"mnn_samples_pca") <- mnn_samples$corrected

# Run batch correction with MNN across Species 
mnn_species <- batchelor::reducedMNN(reducedDim(combined, "PCA"),
                                     batch=as.factor(combined$Species))

reducedDim(combined,"mnn_species_pca") <- mnn_species$corrected

combined
# class: SingleCellExperiment 
# dim: 14391 187974 
# metadata(0):
#     assays(4): merged counts logcounts binomial_deviance_residuals
# rownames(14391): SAMD11 NOC2L ... EIF1AY RPS4Y2
# rowData names(4): gene_name gene_type Symbol.uniq gene_biotype
# colnames(187974): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# colData names(16): batch Sample ... sizeFactor Species
# reducedDimNames(6): PCA TSNE ... mnn_samples mnn_species
# mainExpName: NULL
# altExpNames(0):



# ============= TSNE and UMAP on MNN corrected data ===============

# calculate UMAP and TSNE using GLM-PCA_approx
set.seed(1234)
# combined <- runUMAP(combined,
#                     dimred = "mnn_species",
#                     name = "UMAP_mnn_species",
#                     #BPPARAM = BiocParallel::MulticoreParam(workers=20)
#                     )

combined <- runTSNE(combined,
                    dimred = "mnn_species_pca",
                    name = "TSNE_mnn_species_pca"
                    #BPPARAM = BiocParallel::MulticoreParam(workers=20)
)

pdf(here(plot_dir, "TSNE_mnn_corrected_species_pca.pdf"))
plotReducedDim(combined, dim="TSNE_mnn_species_pca", colour_by = "Species", point_alpha = 0.2)
dev.off()



combined <- runTSNE(combined,
                    dimred = "mnn_samples_pca",
                    name = "TSNE_mnn_samples_pca"
                    #BPPARAM = BiocParallel::MulticoreParam(workers=20)
)

pdf(here(plot_dir, "TSNE_mnn_corrected_samples_pca.pdf"))
plotReducedDim(combined, dim="TSNE_mnn_samples_pca", colour_by = "Species", point_alpha = 0.2)
dev.off()




# ========= Batch correction by samples and species ==============


# Run batch correction with MNN across Samples
mnn_samples <- batchelor::reducedMNN(reducedDim(combined, "PCA"),
                                     batch=as.factor(combined$Sample))

reducedDim(combined,"mnn_samples") <- mnn_samples$corrected

# Run batch correction with MNN across Species 
mnn_species <- batchelor::reducedMNN(reducedDim(combined, "mnn_samples"),
                                     batch=as.factor(combined$Species))

reducedDim(combined,"mnn_samples_species") <- mnn_species$corrected

combined
# class: SingleCellExperiment 
# dim: 14391 187974 
# metadata(0):
#     assays(4): merged counts logcounts binomial_deviance_residuals
# rownames(14391): SAMD11 NOC2L ... EIF1AY RPS4Y2
# rowData names(4): gene_name gene_type Symbol.uniq gene_biotype
# colnames(187974): 1_AAACCCAAGCTAAATG-1 1_AAACCCACAGGTCCCA-1 ...
# 7_TTTGTTGTCGTGGGAA-1 7_TTTGTTGTCTTAATCC-1
# colData names(16): batch Sample ... sizeFactor Species
# reducedDimNames(6): PCA TSNE ... mnn_samples mnn_species
# mainExpName: NULL
# altExpNames(0):



# ============= TSNE and UMAP on MNN corrected data ===============

# calculate UMAP and TSNE using GLM-PCA_approx
set.seed(1234)
# combined <- runUMAP(combined,
#                     dimred = "mnn_species",
#                     name = "UMAP_mnn_species",
#                     #BPPARAM = BiocParallel::MulticoreParam(workers=20)
#                     )

combined <- runTSNE(combined,
                    dimred = "mnn_samples_species",
                    name = "TSNE_mnn_samples_species",
                    BPPARAM = BiocParallel::MulticoreParam(workers=20)
)

pdf(here(plot_dir, "TSNE_mnn_samples_and_species.pdf"))
plotReducedDim(combined, dim="TSNE_mnn_samples_species", colour_by = "Species", point_alpha = 0.2)
dev.off()

