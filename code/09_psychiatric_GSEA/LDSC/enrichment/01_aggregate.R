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
    file = here("processed-data", "09_psychiatric_GSEA", "LDSC","enrichment", "human_DE_results.Rdata")
)

# create a matrix where rows are genes, columns are clusters (domains) and elements are some statistic from DE analysis (eg t-stats)
# we will use the t-statistic for now

enrichment_results <- modeling_results[["enrichment"]]

which(enrichment_results$gene == "HSPA14")
#[1] 9656 9657
enrichment_results <- enrichment_results[-9656,]

# find repeated gene names in enrichment_results$gene
duplicated_genes <- enrichment_results[duplicated(enrichment_results$gene),]$gene

# remove one of the duplicated genes at random for each gene
for (gene in duplicated_genes) {
    enrichment_results <- enrichment_results[-sample(which(enrichment_results$gene == gene), 1),]
}

k <- unique(sce$fine_celltype)


aggregated <- enrichment_results[,c(1:(length(k)-3))]

#remove "t_stat_" prefix from column names
colnames(aggregated) <- gsub("t_stat_", "", colnames(aggregated))

rownames(aggregated) <- enrichment_results$gene

write.table(aggregated, here::here("processed-data", "09_psychiatric_GSEA","LDSC","enrichment", "human_aggregated_de.tsv"), na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
