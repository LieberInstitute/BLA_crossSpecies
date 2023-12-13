install.packages("rafalib")

library("SingleCellExperiment")
library("scuttle")
library("scran")
library("scater")
library("scDblFinder")
#library("jaffelab")
library("batchelor")
library("tidyverse")
library("here")
library("sessioninfo")
library("HDF5Array")
library("rafalib")
library("dplyr")


## Load raw data
load(here("processed-data", "02_build_sce", "sce_baboon_raw_ncbi.rda"), verbose = TRUE)
sce

pd <- colData(sce) %>% as.data.frame()

## save directories
plot_dir = here("plots", "03_quality_control","Baboon","PerCellQC")
processed_dir = here("processed-data","03_quality_control","Baboon","PerCellQC")

##get droplet score filepaths
droplet_paths <- list.files(here("processed-data","03_quality_control","Baboon","droplet_scores"),
                            full.names = TRUE
)

names(droplet_paths) <- gsub("st", "s", gsub("droplet_scores_|.Rdata", "", basename(droplet_paths)))

e.out <- lapply(droplet_paths, function(x) get(load(x)))


# Read the log file
log_file_path <- here('code','03_quality_control','logs','emptyDrops_baboon.txt')
log_text <- readLines(log_file_path, warn = FALSE)
log_text

# Find the lines with 'knee_lower' values
knee_lower_lines <- grep('knee_lower', log_text, value = TRUE)
knee_lower_lines

# Extract the 'knee_lower' values
knee_lower_values <- str_extract(knee_lower_lines, "(?<=knee_lower =)\\d+")
knee_lower_values

knee_lower <- as.numeric(knee_lower_values)
#names(knee_lower) <- c('BR2327','Br8692','Br9021','Br8337', 'Br5273')
knee_lower

## This is weird. doesn't work with my logs. Easier to just merge these scripts? - MT
#knee_lower <- map_dbl(logs, ~ parse_number(.x[grepl("knee_lower =", .x)]))
#names(knee_lower) <- gsub("st", "s", map_chr(logs, ~str_sub(.x[grepl("Running Sample: ", .x)], " ", 2)))

names(knee_lower) <- c(
    "SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1",
    "SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1",
    "SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1",
    "SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1",   
    "SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1",
    "SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1",          
    "SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1"
)
knee_lower


e.out <- e.out[names(knee_lower)]
##n empty table

FDR_cutoff <- 0.001

drop_summary <- stack(map_int(e.out, nrow)) %>%
    rename(total_n = values) %>%
    left_join(stack(map_int(e.out, ~ sum(.x$FDR < FDR_cutoff, na.rm = TRUE))) %>%
                  rename(non_empty = values)) %>%
    select(Sample = ind, total_n, non_empty) %>%
    left_join(stack(knee_lower) %>% rename(Sample = ind, lower_cutoff = values))

write_csv(drop_summary, file = here(processed_dir, "drop_summary.csv"))

drop_summary %>%
    arrange(non_empty)

summary(drop_summary$non_empty)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4867    5028    5190    5190    5351    5512 

#make Louise-style barplot
drop_barplot <- drop_summary %>%
    mutate(empty = total_n - non_empty) %>%
    select(-total_n) %>%
    pivot_longer(!Sample, names_to = "drop_type", values_to = "n_drop") %>%
    ggplot(aes(x = Sample, y = n_drop, fill = drop_type)) +
    geom_col() +
    scale_y_continuous(trans = "log10") +
    #my_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(drop_barplot, filename = here(plot_dir, "drop_barplot.png"), width = 9)


#drop_v_cutoff <- drop_summary %>%
#    left_join(round_info) %>%
#    ggplot(aes(x = lower_cutoff, y = non_empty, color = round)) +
#    geom_point() +
#    my_theme

#ggsave(drop_v_cutoff, filename = here(plot_dir, "drop_v_cutoff.png"))

## Check empty droplet results
map(e.out, ~ addmargins(table(Signif = .x$FDR <= FDR_cutoff, Limited = .x$Limited, useNA = "ifany")))

# $SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   26809       0       0   26809
# TRUE      103    6386       0    6489
# <NA>        0       0 1569037 1569037
# Sum     26912    6386 1569037 1602335
# 
# $SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   32904       0       0   32904
# TRUE      124    6164       0    6288
# <NA>        0       0 1533919 1533919
# Sum     33028    6164 1533919 1573111
# 
# $SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE    2892       0       0    2892
# TRUE     3844    4100       0    7944
# <NA>        0       0 1683928 1683928
# Sum      6736    4100 1683928 1694764
# 
# $SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE     558       0       0     558
# TRUE       61    6170       0    6231
# <NA>        0       0 1701216 1701216
# Sum       619    6170 1701216 1708005
# 
# $SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE     176       0       0     176
# TRUE      220    7780       0    8000
# <NA>        0       0 2219738 2219738
# Sum       396    7780 2219738 2227914
# 
# $SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE     790       0       0     790
# TRUE      160   10186       0   10346
# <NA>        0       0 2534554 2534554
# Sum       950   10186 2534554 2545690
# 
# $SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE     991       0       0     991
# TRUE      841   12223       0   13064
# <NA>        0       0 2573192 2573192
# Sum      1832   12223 2573192 2587247


#### Eliminate empty droplets ####
#e.out.all <- do.call("rbind", e.out)[colnames(sce), ]
e.out.all <- do.call("rbind", e.out)

sce.dropped <- sce[, which(e.out.all$FDR <= 0.001)]

dim(sce.dropped)
<<<<<<< HEAD
# [1] 21091 58362
=======
# 21091 58362
>>>>>>> 46be773 (added in 3 missing mito genes. rebuilt SCE and reran PerCellQC)

# use grep to find genes with "MT-" in name
is.mito <- grep("MT-", rownames(sce.dropped))

sce.dropped <- scuttle::addPerCellQC(
    sce.dropped,
    subsets = list(Mito = is.mito),
    # BPPARAM = BiocParallel::MulticoreParam(4)
)

colnames(colData(sce.dropped))
# [1] "sample_id"             "Sample"                "subject"              
# [4] "species"               "region"                "subregion"            
# [7] "dv_axis"               "Barcode"               "key"                  
# [10] "sum"                   "detected"              "subsets_Mito_sum"     
# [13] "subsets_Mito_detected" "subsets_Mito_percent"  "total" 



# #save drops removed sce
save(sce.dropped,file=here(processed_dir,"sce_drops_removed_baboon.rda"))

load(here(processed_dir, "sce_drops_removed_baboon.rda"), verbose = TRUE)

sce <- sce.dropped

#### Check for low quality nuc ####
# Sample 4 is bimodal and so MAD does not work for it. subsetting to all other samples
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, 
                           nmads = 3, 
                           type = "higher", 
                           batch = sce$Sample,
                           subset=sce$Sample %in% c("SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1", 
                                                    "SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1", 
                                                    "SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1",
                                                    "SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1",
                                                    "SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1",
                                                    "SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1"
                                                    )
                           )

table(sce$high_mito)

# before subset
# FALSE  TRUE 
# 46564 11798 

# after subset
# FALSE  TRUE 
# 42570 15792 

# after fixing mito genes
# FALSE  TRUE 
# 49373  8989 

table(sce$high_mito, sce$Sample)

# SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1
# FALSE                                        6234
# TRUE                                          255
# 
# SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1
# FALSE                                       5940
# TRUE                                         348
# 
# SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1
# FALSE                                       7433
# TRUE                                         511
# 
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1
# FALSE                                       2063
# TRUE                                        4168
# 
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1
# FALSE                                                 6523
# TRUE                                                  1477
# 
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1
# FALSE                                       9850
# TRUE                                         496
# 
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1
# FALSE                                                11330
# TRUE                                                  1734

## low library size
# samples 1 and 3 are bimodal here. subsetting to all other samples
sce$low_lib <- isOutlier(sce$sum, 
                         log = TRUE, 
                         type = "lower", 
                         batch = sce$Sample,
                         subset=sce$Sample %in% c("SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1",
                                                  "SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1",
                                                  "SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1",
                                                  "SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1",
                                                  "SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1"
                                                  )
                         )
table(sce$low_lib)
# before subset
# FALSE  TRUE 
# 57863   499 

# after subset
# FALSE  TRUE 
# 55539  2823 

# after fixed genes
# FALSE  TRUE 
# 55539  2823 

## low detected features
# samples 1 and 3 are also bimodel here. subsettings to all others
sce$low_genes <- isOutlier(sce$detected, 
                           log = TRUE, 
                           type = "lower", 
                           batch = sce$Sample,
                           subset=sce$Sample %in% c("SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1",
                                                   "SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1",
                                                   "SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1",
                                                   "SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1",
                                                   "SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1"
                                                   )
                           )
table(sce$low_genes)
# before subset
# FALSE  TRUE 
# 57494   868  

# after subset
# FALSE  TRUE 
# 53487  4875 

# after fixed genes
# FALSE  TRUE 
# 53487  4875 



## All low sum are also low detected
table(sce$low_lib, sce$low_genes)
#       FALSE  TRUE
# FALSE 53487  2052
# TRUE      0  2823

## Annotate nuc to drop
sce$discard_auto <- sce$high_mito | sce$low_lib | sce$low_genes

table(sce$discard_auto)
# FALSE  TRUE 
# 47039 11323 



qc_t <- addmargins(table(sce$Sample, sce$discard_auto))
write_csv(data.frame(qc_t), file = here(processed_dir, "discarded_summary.csv"))

qc_t
#                                                      FALSE  TRUE   Sum
# FALSE  TRUE   Sum
# SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1           6166   323  6489
# SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1            5802   486  6288
# SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1            7341   603  7944
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1            1922  4309  6231
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1  5596  2404  8000
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1            9763   583 10346
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1 10449  2615 13064
# Sum                                                  47039 11323 58362



round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)

# SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1           95.0   5.0 100.0
# SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1            92.3   7.7 100.0
# SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1            92.4   7.6 100.0
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1            30.8  69.2 100.0
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1  70.0  30.0 100.0
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1            94.4   5.6 100.0
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1  80.0  20.0 100.0
# Sum                                                   80.6  19.4 100.0


#### QC plots ####

# Mito rate
png(here(plot_dir, "violin_mito_percent_afterSubset_fixedGenes.png"), height=7.5, width=7.5, units="in", res=1200)
plotColData(sce, x = "subsets_Mito_percent", y = "Sample", colour_by = "high_mito") +
    ggtitle("Mito Percent") +
    theme(axis.text.x = element_text(angle = 45))
dev.off()

# low sum
png(here(plot_dir, "violin_sum_afterSubset_fixedGenes.png"), height=7.5, width=7.5, units="in", res=1200)
plotColData(sce, x = "sum", y = "Sample", colour_by = "low_lib") +
    scale_y_log10() +
    ggtitle("Total UMIs")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# low detected
png(here(plot_dir, "violin_detected_afterSubset_fixedGenes.png"), height=7.5, width=7.5, units="in", res=1200)
plotColData(sce, x = "detected", y = "Sample", colour_by = "low_genes") +
    scale_y_log10() +
    ggtitle("Detected genes")
dev.off()





# Mito rate vs n detected features
png(here(plot_dir, "scatter_detect_vs_mit.png"), , height=7.5, width=7.5, units="in", res=1200)
plotColData(sce,
            x = "detected", y = "subsets_Mito_percent",
            colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5
) + 
    ggtitle("Mito percent") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# Detected features vs total count
png(here(plot_dir, "scatter_detect_vs_sum.png"), , height=7.5, width=7.5, units="in", res=1200)
plotColData(sce,
            x = "sum", y = "detected",
            colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5
) + 
    ggtitle("Sum detected") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



#### Doublet detection ####
## To speed up, run on sample-level top-HDGs - just take top 2000
set.seed(328)

colData(sce)$doubletScore <- NA

for (i in splitit(sce$Sample)) {
    sce_temp <- sce[, i]
    ## To speed up, run on sample-level top-HVGs - just take top 1000
    normd <- logNormCounts(sce_temp)
    geneVar <- modelGeneVar(normd)
    topHVGs <- getTopHVGs(geneVar, n = 2000)
    
    dbl_dens <- computeDoubletDensity(normd, subset.row = topHVGs)
    colData(sce)$doubletScore[i] <- dbl_dens
}

summary(sce$doubletScore)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0400  0.2000  0.7542  0.6480 21.1025 

quantile(sce$doubletScore, probs=seq(0,1,by=0.01),3)


## Visualize doublet scores ##

dbl_df <- colData(sce) %>%
    as.data.frame() %>%
    dplyr::select(Sample, doubletScore)

dbl_box_plot <- dbl_df %>%
    ggplot(aes(x =reorder(Sample, doubletScore, FUN = median) , y = doubletScore)) +
    geom_boxplot() +
    labs(x = "Sample") +
    geom_hline(yintercept = 2.75, color = "red", linetype = "dashed") +
    coord_flip() +
    theme(text = element_text(size = 30)) 

ggsave(dbl_box_plot, filename = here(plot_dir, "doublet_scores_boxplot.png"), width=20, height=25)

#dbl_density_plot <- dbl_df %>%
#  ggplot(aes(x = doubletScore)) +
#  geom_density() +
#  labs(x = "doublet score") +
#  facet_grid(Sample ~ .) +
#  theme_bw() + 
#  theme(text = element_text(size = 30)) 

#ggsave(dbl_density_plot, filename = here(plot_dir, "doublet_scores_density.png"), height = 17)

dbl_df %>%
    group_by(Sample) %>%
    summarize(
        median = median(doubletScore),
        q95 = quantile(doubletScore, .95),
        drop = sum(doubletScore >= 2.75),
        drop_percent = 100 * drop / n()
    )

#   Sample                                         median   q95  drop drop_percent
# <chr>                                           <dbl> <dbl> <int>        <dbl>
# 1 SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1     0.389  2.44   275         4.24
# 2 SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1      0.365  2.00   210         3.34
# 3 SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1      0.159  3.50   467         5.88
# 4 SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1      0.411  4.01   507         8.14
# 5 SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI…  0.16   3.26   620         7.75
# 6 SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1      0.16   3.88   607         5.87
# 7 SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI…  0.1    4.94  1159         8.87




table(sce$discard_auto, sce$doubletScore >= 2.75)
#       FALSE  TRUE
# FALSE 43715  3324
# TRUE  10802   521


#### Save clean data as HDF5 file  ####
load(here(processed_dir, "sce_no_empty_droplets.Rdata"))

sce <- sce[, !sce$discard_auto]
dim(sce)
# [1] 21091 47039

save(sce,file=here(processed_dir,"sce_post_qc_baboon.rda"))


# sgejobs::job_single('03_droplet_qc', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 03_droplet_qc.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 Patched (2023-11-13 r85524)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-12-13
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3.x/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                * 1.4-5     2016-07-21 [2] CRAN (R 4.3.2)
# batchelor            * 1.18.0    2023-10-24 [2] Bioconductor
# beachmat               2.18.0    2023-10-24 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.3.2)
# Biobase              * 2.62.0    2023-10-24 [2] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [2] Bioconductor
# BiocIO                 1.12.0    2023-10-24 [2] Bioconductor
# BiocNeighbors          1.20.0    2023-10-24 [2] Bioconductor
# BiocParallel           1.36.0    2023-10-24 [2] Bioconductor
# BiocSingular           1.18.0    2023-10-24 [2] Bioconductor
# Biostrings             2.70.1    2023-10-25 [2] Bioconductor
# bit                    4.0.5     2022-11-15 [2] CRAN (R 4.3.2)
# bit64                  4.0.5     2020-08-30 [2] CRAN (R 4.3.2)
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.2)
# bluster                1.12.0    2023-10-24 [2] Bioconductor
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.2)
# cluster                2.1.4     2022-08-22 [3] CRAN (R 4.3.2)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.2)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.2)
# cowplot                1.1.1     2020-12-30 [2] CRAN (R 4.3.2)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.2)
# data.table             1.14.8    2023-02-17 [2] CRAN (R 4.3.2)
# DelayedArray         * 0.28.0    2023-10-24 [2] Bioconductor
# DelayedMatrixStats     1.24.0    2023-10-24 [2] Bioconductor
# dplyr                * 1.1.3     2023-09-03 [2] CRAN (R 4.3.2)
# dqrng                  0.3.1     2023-08-30 [2] CRAN (R 4.3.2)
# DropletUtils         * 1.22.0    2023-10-24 [2] Bioconductor
# edgeR                  4.0.1     2023-10-29 [2] Bioconductor
# fansi                  1.0.5     2023-10-08 [2] CRAN (R 4.3.2)
# farver                 2.1.1     2022-07-06 [2] CRAN (R 4.3.2)
# forcats              * 1.0.0     2023-01-29 [2] CRAN (R 4.3.2)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.2)
# GenomeInfoDb         * 1.38.1    2023-11-08 [2] Bioconductor
# GenomeInfoDbData       1.2.11    2023-11-15 [2] Bioconductor
# GenomicAlignments      1.38.0    2023-10-24 [2] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [2] Bioconductor
# ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.3.2)
# ggplot2              * 3.4.4     2023-10-12 [2] CRAN (R 4.3.2)
# ggrepel                0.9.4     2023-10-13 [2] CRAN (R 4.3.2)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.2)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.3.2)
# gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.2)
# HDF5Array            * 1.30.0    2023-10-24 [2] Bioconductor
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.2)
# hms                    1.1.3     2023-03-21 [2] CRAN (R 4.3.2)
# igraph                 1.5.1     2023-08-10 [2] CRAN (R 4.3.2)
# IRanges              * 2.36.0    2023-10-24 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.2)
# jsonlite               1.8.7     2023-06-29 [2] CRAN (R 4.3.2)
# labeling               0.4.3     2023-08-29 [2] CRAN (R 4.3.2)
# lattice                0.22-5    2023-10-24 [3] CRAN (R 4.3.2)
# lifecycle              1.0.4     2023-11-07 [2] CRAN (R 4.3.2)
# limma                  3.58.1    2023-10-31 [2] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [2] CRAN (R 4.3.2)
# lubridate            * 1.9.3     2023-09-27 [2] CRAN (R 4.3.2)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.2)
# MASS                   7.3-60    2023-05-04 [3] CRAN (R 4.3.2)
# Matrix               * 1.6-4     2023-11-30 [1] CRAN (R 4.3.2)
# MatrixGenerics       * 1.14.0    2023-10-24 [2] Bioconductor
# matrixStats          * 1.1.0     2023-11-07 [2] CRAN (R 4.3.2)
# metapod                1.10.0    2023-10-24 [2] Bioconductor
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.2)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.2)
# purrr                * 1.0.2     2023-08-10 [2] CRAN (R 4.3.2)
# R.methodsS3            1.8.2     2022-06-13 [2] CRAN (R 4.3.2)
# R.oo                   1.25.0    2022-06-12 [2] CRAN (R 4.3.2)
# R.utils                2.12.2    2022-11-11 [2] CRAN (R 4.3.2)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.2)
# rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 4.3.2)
# ragg                   1.2.6     2023-10-10 [2] CRAN (R 4.3.2)
# RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.2)
# RCurl                  1.98-1.13 2023-11-02 [2] CRAN (R 4.3.2)
# readr                * 2.1.4     2023-02-10 [2] CRAN (R 4.3.2)
# ResidualMatrix         1.12.0    2023-10-24 [2] Bioconductor
# restfulr               0.0.15    2022-06-16 [2] CRAN (R 4.3.2)
# rhdf5                * 2.46.0    2023-10-24 [2] Bioconductor
# rhdf5filters           1.14.1    2023-11-06 [2] Bioconductor
# Rhdf5lib               1.24.0    2023-10-24 [2] Bioconductor
# rjson                  0.2.21    2022-01-09 [2] CRAN (R 4.3.2)
# rlang                  1.1.2     2023-11-04 [2] CRAN (R 4.3.2)
# rprojroot              2.0.4     2023-11-05 [2] CRAN (R 4.3.2)
# Rsamtools              2.18.0    2023-10-24 [2] Bioconductor
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.3.2)
# rtracklayer          * 1.62.0    2023-10-24 [2] Bioconductor
# S4Arrays             * 1.2.0     2023-10-24 [2] Bioconductor
# S4Vectors            * 0.40.1    2023-10-26 [2] Bioconductor
# ScaledMatrix           1.10.0    2023-10-24 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.2)
# scater               * 1.30.0    2023-10-24 [2] Bioconductor
# scDblFinder          * 1.16.0    2023-10-24 [1] Bioconductor
# scran                * 1.30.0    2023-10-24 [2] Bioconductor
# scuttle              * 1.12.0    2023-10-24 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.2)
# SingleCellExperiment * 1.24.0    2023-10-24 [2] Bioconductor
# SparseArray          * 1.2.2     2023-11-07 [2] Bioconductor
# sparseMatrixStats      1.14.0    2023-10-24 [2] Bioconductor
# statmod                1.5.0     2023-01-06 [2] CRAN (R 4.3.2)
# stringi                1.8.1     2023-11-13 [2] CRAN (R 4.3.2)
# stringr              * 1.5.1     2023-11-14 [2] CRAN (R 4.3.2)
# SummarizedExperiment * 1.32.0    2023-10-24 [2] Bioconductor
# systemfonts            1.0.5     2023-10-09 [2] CRAN (R 4.3.2)
# textshaping            0.3.7     2023-10-09 [2] CRAN (R 4.3.2)
# tibble               * 3.2.1     2023-03-20 [2] CRAN (R 4.3.2)
# tidyr                * 1.3.0     2023-01-24 [2] CRAN (R 4.3.2)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.2)
# tidyverse            * 2.0.0     2023-02-22 [2] CRAN (R 4.3.2)
# timechange             0.2.0     2023-01-11 [2] CRAN (R 4.3.2)
# tzdb                   0.4.0     2023-05-12 [2] CRAN (R 4.3.2)
# utf8                   1.2.4     2023-10-22 [2] CRAN (R 4.3.2)
# vctrs                  0.6.4     2023-10-12 [2] CRAN (R 4.3.2)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.3.2)
# viridis                0.6.4     2023-07-22 [2] CRAN (R 4.3.2)
# viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.3.2)
# vroom                  1.6.4     2023-10-02 [2] CRAN (R 4.3.2)
# withr                  2.5.2     2023-10-30 [2] CRAN (R 4.3.2)
# xgboost                1.7.5.1   2023-03-30 [2] CRAN (R 4.3.2)
# XML                    3.99-0.15 2023-11-02 [2] CRAN (R 4.3.2)
# XVector                0.42.0    2023-10-24 [2] Bioconductor
# yaml                   2.3.7     2023-01-23 [2] CRAN (R 4.3.2)
# zlibbioc               1.48.0    2023-10-24 [2] Bioconductor
# 
# [1] /users/mtotty/R/4.3.x
# [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# > 