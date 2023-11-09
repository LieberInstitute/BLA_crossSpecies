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
load(here("processed-data", "02_build_sce", "sce_baboon_raw_no_sample_info.rda"), verbose = TRUE)
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

# $`34ac_scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   10567       0       0   10567
# TRUE       49    4818       0    4867
# <NA>        0       0 1359784 1359784
# Sum     10616    4818 1359784 1375218
# 
# $`35ac_scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE    6678       0       0    6678
# TRUE      177    5335       0    5512
# <NA>        0       0 1053461 1053461
# Sum      6855    5335 1053461 1065651
# 
# $`3c-AMYBLA`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   14453       0       0   14453
# TRUE      248    4518       0    4766
# <NA>        0       0 1721248 1721248
# Sum     14701    4518 1721248 1740467
# 
# $`4c-AMYBLA`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   11401       0       0   11401
# TRUE      111    4887       0    4998
# <NA>        0       0 1773286 1773286
# Sum     11512    4887 1773286 1789685
# 
# $`5c-AMYBLA`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   15192       0       0   15192
# TRUE       83    5176       0    5259
# <NA>        0       0 1840574 1840574
# Sum     15275    5176 1840574 1861025


#### Eliminate empty droplets ####
#e.out.all <- do.call("rbind", e.out)[colnames(sce), ]
e.out.all <- do.call("rbind", e.out)

sce.dropped <- sce[, which(e.out.all$FDR <= 0.001)]

dim(sce.dropped)
# [1]  21369 134685


#### Compute QC metrics ####
sce.dropped <- scuttle::addPerCellQC(
    sce.dropped,
    subsets = list(Mito = which(seqnames(sce) == "MT")),
    # BPPARAM = BiocParallel::MulticoreParam(4)
)

colnames(colData(sce.dropped))
# [1] "sample_id"             "Sample"                "subject"              
# [4] "species"               "region"                "subregion"            
# [7] "dv_axis"               "Barcode"               "key"                  
# [10] "sum"                   "detected"              "subsets_Mito_sum"     
# [13] "subsets_Mito_detected" "subsets_Mito_percent"  "total" 



# #save drops removed sce
save(sce.dropped,file=here(processed_dir,"sce_drops_removed.rda"))

load(here(processed_dir, "sce_drops_removed.rda"), verbose = TRUE)

sce <- sce.dropped

#### Check for low quality nuc ####
## High mito
# sce$high.mito.sample ## standard name?
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$Sample)

table(sce$high_mito)
# FALSE   TRUE 
# 109967  24718 

table(sce$high_mito, sce$Sample)


# LIB210527RC_AB_1A LIB210527RC_AB_1B
# FALSE              4184              3856
# TRUE                799               752
# 
# VC_snRNAseq_12_Animal2_Central_Nucleus
# FALSE                                   3434
# TRUE                                     674
# 
# VC_snRNAseq_12_Animal3_Central_Nucleus
# FALSE                                   4548
# TRUE                                     974
# 
# VC_snRNAseq_12_Animal4_Central_Nucleus
# FALSE                                   2130
# TRUE                                    2048
# 
# VC_snRNAseq_12_Animal5_Accessory_Basal__AP2_AP3_
# FALSE                                             2108
# TRUE                                               278
# 
# VC_snRNAseq_12_Animal5_Basal__AP1_ VC_snRNAseq_12_Animal5_Basal__AP2_
# FALSE                               2575                               3280
# TRUE                                 340                                665
# 
# VC_snRNAseq_12_Animal5_Lateral__AP1_
# FALSE                                  470
# TRUE                                    90
# 
# VC_snRNAseq_12_Animal5_Lateral__AP2_
# FALSE                                 2855
# TRUE                                   395
# 
# VC_snRNAseq_13_Animal5_Central_Nucleus__AP3_ VC_snRNAseq_8_Animal3_AB__
# FALSE                                         1825                       2581
# TRUE                                           183                        298
# 
# VC_snRNAseq_8_Animal3_BD VC_snRNAseq_8_Animal3_Bd_Bv
# FALSE                     2147                        3507
# TRUE                       360                         633
# 
# VC_snRNAseq_8_Animal3_BV1 VC_snRNAseq_8_Animal3_BV2
# FALSE                      3912                      3670
# TRUE                        255                       266
# 
# VC_snRNAseq_8_Animal3_LD VC_snRNAseq_8_Animal3_LV
# FALSE                     3773                     2350
# TRUE                       558                      260
# 
# VC_snRNAseq_8_Animal3_LV3 VC_snRNAseq_9_Animal4_AB
# FALSE                      1788                     3591
# TRUE                        443                      616
# 
# VC_snRNAseq_9_Animal4_B-Comb VC_snRNAseq_9_Animal4_BD
# FALSE                         3599                     4073
# TRUE                           799                     1313
# 
# VC_snRNAseq_9_Animal4_BV VC_snRNAseq_9_Animal4_L-Comb
# FALSE                     3563                         3834
# TRUE                      1401                         1366
# 
# VC_snRNAseq_9_Animal4_LD VC_snRNAseq_9_Animal4_LV1
# FALSE                     3819                      2216
# TRUE                      2863                      1779
# 
# VC_snRNAseq_9_Animal4_LV2 VC_snRNAseq-7_AccBasalAP1AP2_-7
# FALSE                      2702                            3446
# TRUE                        487                             930
# 
# VC_snRNAseq-7_BasalDorsalAP1_-3A VC_snRNAseq-7_BasalDorsalAP1_-3B
# FALSE                             3865                             3497
# TRUE                               439                              464
# 
# VC_snRNAseq-7_BasalVentralAP1_-4 VC_snRNAseq-7_BasalVentralAP2_-5
# FALSE                             1962                             3443
# TRUE                               245                              568
# 
# VC_snRNAseq-7_LateralDorsalAP1AP2_-8 VC_snRNAseq-7_LateralVentralAP1_-1
# FALSE                                 4276                               2764
# TRUE                                   468                                316
# 
# VC_snRNAseq-7_LateralVentralAP2_-2
# FALSE                               4324
# TRUE                                 393

## low library size
sce$low_lib <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sample)
table(sce$low_lib)
# FALSE  TRUE 
# 57886   476 

## low detected features
# sce$qc.detected
sce$low_genes <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$Sample)
table(sce$low_genes)
# FALSE  TRUE 
# 57564   798 



## All low sum are also low detected
table(sce$low_lib, sce$low_genes)
# FALSE  TRUE
# FALSE 57564   322
# TRUE      0   476


## Annotate nuc to drop
sce$discard_auto <- sce$high_mito | sce$low_lib | sce$low_genes

table(sce$discard_auto)
# FALSE   TRUE 
# 109142  25543 



qc_t <- addmargins(table(sce$Sample, sce$discard_auto))
write_csv(data.frame(qc_t), file = here(processed_dir, "discarded_summary.csv"))

qc_t
#                                                   FALSE   TRUE    Sum
# FALSE  TRUE   Sum
# SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1           6389   100  6489
# SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1            6091   197  6288
# SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1            7796   148  7944
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1            6061   170  6231
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1  8000     0  8000
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1           10163   183 10346
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1 13064     0 13064
# Sum                                                  57564   798 58362



round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)

# FALSE  TRUE   Sum
# SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1           98.5   1.5 100.0
# SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1            96.9   3.1 100.0
# SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1            98.1   1.9 100.0
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1            97.3   2.7 100.0
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1 100.0   0.0 100.0
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1            98.2   1.8 100.0
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1 100.0   0.0 100.0
# Sum                                                   98.6   1.4 100.0



#### QC plots ####
pdf(height=15, width=7.5, here(plot_dir, "QC_violin_plots.pdf"))
## Mito rate
plotColData(sce, x = "subsets_Mito_percent", y = "Sample", colour_by = "high_mito") +
    ggtitle("Mito Percent") +
    theme(axis.text.x = element_text(angle = 45))

# ## low sum
plotColData(sce, x = "sum", y = "Sample", colour_by = "low_lib") +
    scale_y_log10() +
    ggtitle("Total UMIs")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ## low detected
plotColData(sce, x = "detected", y = "Sample", colour_by = "low_genes") +
    scale_y_log10() +
    ggtitle("Detected genes")

dev.off()





# Mito rate vs n detected features
pdf(height=7.5, width=7.5, here(plot_dir, "QC_scatter_plots.pdf"))
plotColData(sce,
            x = "detected", y = "subsets_Mito_percent",
            colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5
) + 
    ggtitle("Mito percent") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Detected features vs total count
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
# 0.0000  0.1168  0.2840  0.5284  0.5815 21.6952 

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

# Sample    median   q95  drop drop_percent
# <chr>      <dbl> <dbl> <int>        <dbl>
# 1 34ac_scp   0.224  1.09    90         1.85
# 2 35ac_scp   0.287  1.27   120         2.18
# 3 3c-AMYBLA  0.429  1.61    59         1.24
# 4 4c-AMYBLA  0.250  1.24   112         2.24
# 5 5c-AMYBLA  0.316  1.53    93         1.77




table(sce$discard_auto, sce$doubletScore >= 2.75)
#       FALSE  TRUE
# FALSE 20821   447
# TRUE   4107    27


#### Save clean data as HDF5 file  ####
load(here(processed_dir, "sce_no_empty_droplets.Rdata"))

sce <- sce[, !sce$discard_auto]
dim(sce)
# [1] 36601 21268

save(sce,file=here(processed_dir,"sce_post_qc.rda"))


# sgejobs::job_single('03_droplet_qc', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 03_droplet_qc.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
