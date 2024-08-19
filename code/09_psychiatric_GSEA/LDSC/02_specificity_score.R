# originally copied from https://github.com/LieberInstitute/spatialdACC/blob/main/code/17_LDSC/spatial/02_specificity_score.R

setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

# Load libraries
library(here)

is_top_10_percent <- function(column) {
    top_10_percent_threshold <- quantile(column, 0.9)
    as.numeric(column >= top_10_percent_threshold)
}

dat <- read.table(here::here("processed-data", "17_LDSC", "visium_aggregated_cpm.tsv"),header=T)

dat.norm <- t(apply(dat, 1, function(x){x/sum(x)}))
res <- apply(dat.norm, 2, is_top_10_percent)
rownames(res) <- rownames(dat.norm)
write.csv(res,here::here("processed-data", "17_LDSC", "visium_score.csv"))