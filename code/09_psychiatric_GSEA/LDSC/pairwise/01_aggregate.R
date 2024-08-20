# originally copied from https://github.com/LieberInstitute/spatialdACC/blob/main/code/17_LDSC/spatial/01_aggregate.R

# Load libraries
library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)
library(here)
library(scuttle)
library(edgeR)

sce <- readRDS(here("processed-data", "09_psychiatric_GSEA", "sce.human_all_genes.rds"))

modeling_results <- registration_wrapper(
    sce,
    var_registration = "fine_celltype",
    var_sample_id = "Sample",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

###save modeling results list
save(
    modeling_results,
    file = here("processed-data", "09_psychiatric_GSEA", "LDSC","pairwise", "human_DE_results.Rdata")
)

# create a matrix where rows are genes, columns are clusters (domains) and elements are some statistic from DE analysis (eg t-stats)
# we will use the t-statistic for now

pairwise_results <- modeling_results[["pairwise"]]

which(pairwise_results$gene == "HSPA14")
#[1] 9656 9657
pairwise_results <- pairwise_results[-9656,]

# find repeated gene names in pairwise_results$gene
duplicated_genes <- pairwise_results[duplicated(pairwise_results$gene),]$gene

# remove one of the duplicated genes at random for each gene
for (gene in duplicated_genes) {
    pairwise_results <- pairwise_results[-sample(which(pairwise_results$gene == gene), 1),]
}

k <- unique(sce$fine_celltype)


aggregated <- pairwise_results[,c(1:length(k))]

#remove "t_stat_" prefix from column names
colnames(aggregated) <- gsub("t_stat_", "", colnames(aggregated))

rownames(aggregated) <- pairwise_results$gene

write.table(aggregated, here::here("processed-data", "09_psychiatric_GSEA","LDSC","pairwise", "human_aggregated_de.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
