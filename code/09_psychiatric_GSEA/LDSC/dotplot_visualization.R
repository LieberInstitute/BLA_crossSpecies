# originally copied from https://github.com/LieberInstitute/spatialdACC/blob/main/code/17_LDSC/spatial/dotplot_visualization.R

setwd('/dcs04/lieber/marmaypag/spatialdACC_LIBD4125/spatialdACC/')

library(ggplot2)
library(RColorBrewer)
library(here)

###load LDSC results
ldsc_results <- read.csv(file=here::here('code','17_LDSC',
                         'spatial','ldsc_results.csv'))

##########dotplots#############
###make -log10FDR column
ldsc_results$log10fdr <- -log10(ldsc_results$FDR)

####to make nmf only###
#ldsc_results <- ldsc_results[ldsc_results$cell %in% paste0('nmf',c(1:100)),]

####to remove FDR > 0.05
#ldsc_results<-ldsc_results[ldsc_results$FDR<0.05,]

###plot
pdf(here('plots','17_LDSC','spatial','ldsc_results.pdf'),width=10,height=10)

ggplot(ldsc_results, aes(x = cell, y = trait, size = log10fdr, color = Coefficient_z.score)) +
    geom_point() +
    scale_size_continuous(range = c(0, 10)) + # Set minimum size to zero for the size scale
    scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0,
                          name = "Coefficient\n(z-score)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size = "-log10(FDR)",x='group')

ggplot(ldsc_results, aes(x = cell, y = trait, size = ifelse(FDR > 0.05, NA, log10fdr), color = Coefficient_z.score)) +
    geom_point() +
    scale_size_continuous(range = c(0, 10)) + # Set minimum size to zero for the size scale
    scale_color_gradient2(low = 'blue', high = 'red', midpoint = 0,
                          name = "Coefficient\n(z-score)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(size = "-log10(FDR)",x='group',
         caption = "removed FDR > 0.05")

dev.off()