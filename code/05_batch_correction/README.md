# BLA cross-species project: 05_batch_correction

The scripts in this directory perform 1) batch correction across samples and 2) across species using Seurat. This may be somewhat redundant but I found (visually) that this provided better cluster mixing. 

<batch_correction_seurat_samples.R> should be ran prior to <batch_correction_seurat_species>. 

The <comparing_method> scripts conducted initial batch correction comparisons between Seurat, MNN, and Harmony. I've kept them here for archival purposes. These scripts are not essential to the reproduction of the final results. 