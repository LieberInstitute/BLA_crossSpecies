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
load(here("processed-data", "02_build_sce", "sce_macaque_raw_no_sample_info.rda"), verbose = TRUE)
sce

pd <- colData(sce) %>% as.data.frame()

## save directories
## save directories
plot_dir = here("plots", "03_quality_control","Macaque","PerCellQC")
processed_dir = here("processed-data","03_quality_control","Macaque","PerCellQC")

##get droplet score filepaths
droplet_paths <- list.files(here("processed-data","03_quality_control","Macaque","droplet_scores"),
                            full.names = TRUE
)

names(droplet_paths) <- gsub("st", "s", gsub("droplet_scores_|.Rdata", "", basename(droplet_paths)))

e.out <- lapply(droplet_paths, function(x) get(load(x)))


# Read the log file
log_file_path <- here('code','03_quality_control','logs','emptyDrops_macaque.txt')
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
  "VC_snRNAseq-7_LateralVentralAP1_-1",
  "VC_snRNAseq-7_LateralVentralAP2_-2",
  "VC_snRNAseq-7_BasalDorsalAP1_-3A",
  "VC_snRNAseq-7_BasalDorsalAP1_-3B",
  "VC_snRNAseq-7_BasalVentralAP1_-4",
  "VC_snRNAseq-7_BasalVentralAP2_-5",
  "VC_snRNAseq-7_AccBasalAP1AP2_-7",
  "VC_snRNAseq-7_LateralDorsalAP1AP2_-8",
  "VC_snRNAseq_9_Animal4_LV2",
  "VC_snRNAseq_9_Animal4_LV1",
  "VC_snRNAseq_9_Animal4_LD",
  "VC_snRNAseq_9_Animal4_L-Comb",
  "VC_snRNAseq_9_Animal4_BV",
  "VC_snRNAseq_9_Animal4_BD",
  "VC_snRNAseq_9_Animal4_B-Comb",
  "VC_snRNAseq_9_Animal4_AB",
  "VC_snRNAseq_8_Animal3_LV",
  "VC_snRNAseq_8_Animal3_LV3",
  "VC_snRNAseq_8_Animal3_LD",
  "VC_snRNAseq_8_Animal3_BV2",
  "VC_snRNAseq_8_Animal3_BV1",
  "VC_snRNAseq_8_Animal3_BD",
  "VC_snRNAseq_8_Animal3_AB__",
  "VC_snRNAseq_8_Animal3_Bd_Bv",
  "LIB210527RC_AB_1A",
  "LIB210527RC_AB_1B",
  "VC_snRNAseq_12_Animal2_Central_Nucleus",
  "VC_snRNAseq_12_Animal3_Central_Nucleus",
  "VC_snRNAseq_12_Animal4_Central_Nucleus",
  "VC_snRNAseq_12_Animal5_Accessory_Basal__AP2_AP3_",
  "VC_snRNAseq_12_Animal5_Basal__AP1_",
  "VC_snRNAseq_12_Animal5_Basal__AP2_",
  "VC_snRNAseq_12_Animal5_Lateral__AP1_",
  "VC_snRNAseq_12_Animal5_Lateral__AP2_",
  "VC_snRNAseq_13_Animal5_Central_Nucleus__AP3_"
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

which(e.out$FDR <= 0.001)

#### Eliminate empty droplets ####
e.out.all <- do.call("rbind", e.out)[colnames(sce), ]
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
# FALSE   TRUE 
# 134514    171 

## low detected features
# sce$qc.detected
sce$low_genes <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$Sample)
table(sce$low_genes)
# FALSE   TRUE 
# 133322   1363 



## All low sum are also low detected
table(sce$low_lib, sce$low_genes)
#        FALSE   TRUE
# FALSE 133322   1192
# TRUE       0    171


## Annotate nuc to drop
sce$discard_auto <- sce$high_mito | sce$low_lib | sce$low_genes

table(sce$discard_auto)
# FALSE   TRUE 
# 109142  25543 


# ======= Additionally exclude base on minimum thresholds ========

qc.lib <- sce$sum < 600
qc.genes <- sce$detected < 500

sce$discard_minimum <- qc.lib | qc.genes 

table(sce$discard_minimum)


# combine auto and minimum discards
sce$discard <- sce$discard_auto | sce$discard_minimum

table(sce$discard)



qc_t <- addmargins(table(sce$Sample, sce$discard_auto))
write_csv(data.frame(qc_t), file = here(processed_dir, "discarded_summary.csv"))

qc_t
#                                                   FALSE   TRUE    Sum
# LIB210527RC_AB_1A                                  4184    799   4983
# LIB210527RC_AB_1B                                  3856    752   4608
# VC_snRNAseq_12_Animal2_Central_Nucleus             3426    682   4108
# VC_snRNAseq_12_Animal3_Central_Nucleus             4500   1022   5522
# VC_snRNAseq_12_Animal4_Central_Nucleus             2128   2050   4178
# VC_snRNAseq_12_Animal5_Accessory_Basal__AP2_AP3_   2097    289   2386
# VC_snRNAseq_12_Animal5_Basal__AP1_                 2575    340   2915
# VC_snRNAseq_12_Animal5_Basal__AP2_                 3264    681   3945
# VC_snRNAseq_12_Animal5_Lateral__AP1_                457    103    560
# VC_snRNAseq_12_Animal5_Lateral__AP2_               2832    418   3250
# VC_snRNAseq_13_Animal5_Central_Nucleus__AP3_       1802    206   2008
# VC_snRNAseq_8_Animal3_AB__                         2581    298   2879
# VC_snRNAseq_8_Animal3_BD                           2147    360   2507
# VC_snRNAseq_8_Animal3_Bd_Bv                        3507    633   4140
# VC_snRNAseq_8_Animal3_BV1                          3912    255   4167
# VC_snRNAseq_8_Animal3_BV2                          3670    266   3936
# VC_snRNAseq_8_Animal3_LD                           3773    558   4331
# VC_snRNAseq_8_Animal3_LV                           2343    267   2610
# VC_snRNAseq_8_Animal3_LV3                          1788    443   2231
# VC_snRNAseq_9_Animal4_AB                           3343    864   4207
# VC_snRNAseq_9_Animal4_B-Comb                       3599    799   4398
# VC_snRNAseq_9_Animal4_BD                           4073   1313   5386
# VC_snRNAseq_9_Animal4_BV                           3563   1401   4964
# VC_snRNAseq_9_Animal4_L-Comb                       3834   1366   5200
# VC_snRNAseq_9_Animal4_LD                           3793   2889   6682
# VC_snRNAseq_9_Animal4_LV1                          2133   1862   3995
# VC_snRNAseq_9_Animal4_LV2                          2450    739   3189
# VC_snRNAseq-7_AccBasalAP1AP2_-7                    3440    936   4376
# VC_snRNAseq-7_BasalDorsalAP1_-3A                   3865    439   4304
# VC_snRNAseq-7_BasalDorsalAP1_-3B                   3497    464   3961
# VC_snRNAseq-7_BasalVentralAP1_-4                   1940    267   2207
# VC_snRNAseq-7_BasalVentralAP2_-5                   3420    591   4011
# VC_snRNAseq-7_LateralDorsalAP1AP2_-8               4276    468   4744
# VC_snRNAseq-7_LateralVentralAP1_-1                 2751    329   3080
# VC_snRNAseq-7_LateralVentralAP2_-2                 4323    394   4717
# Sum                                              109142  25543 134685



round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)

#                                                   FALSE  TRUE   Sum
# LIB210527RC_AB_1A                                 84.0  16.0 100.0
# LIB210527RC_AB_1B                                 83.7  16.3 100.0
# VC_snRNAseq_12_Animal2_Central_Nucleus            83.4  16.6 100.0
# VC_snRNAseq_12_Animal3_Central_Nucleus            81.5  18.5 100.0
# VC_snRNAseq_12_Animal4_Central_Nucleus            50.9  49.1 100.0
# VC_snRNAseq_12_Animal5_Accessory_Basal__AP2_AP3_  87.9  12.1 100.0
# VC_snRNAseq_12_Animal5_Basal__AP1_                88.3  11.7 100.0
# VC_snRNAseq_12_Animal5_Basal__AP2_                82.7  17.3 100.0
# VC_snRNAseq_12_Animal5_Lateral__AP1_              81.6  18.4 100.0
# VC_snRNAseq_12_Animal5_Lateral__AP2_              87.1  12.9 100.0
# VC_snRNAseq_13_Animal5_Central_Nucleus__AP3_      89.7  10.3 100.0
# VC_snRNAseq_8_Animal3_AB__                        89.6  10.4 100.0
# VC_snRNAseq_8_Animal3_BD                          85.6  14.4 100.0
# VC_snRNAseq_8_Animal3_Bd_Bv                       84.7  15.3 100.0
# VC_snRNAseq_8_Animal3_BV1                         93.9   6.1 100.0
# VC_snRNAseq_8_Animal3_BV2                         93.2   6.8 100.0
# VC_snRNAseq_8_Animal3_LD                          87.1  12.9 100.0
# VC_snRNAseq_8_Animal3_LV                          89.8  10.2 100.0
# VC_snRNAseq_8_Animal3_LV3                         80.1  19.9 100.0
# VC_snRNAseq_9_Animal4_AB                          79.5  20.5 100.0
# VC_snRNAseq_9_Animal4_B-Comb                      81.8  18.2 100.0
# VC_snRNAseq_9_Animal4_BD                          75.6  24.4 100.0
# VC_snRNAseq_9_Animal4_BV                          71.8  28.2 100.0
# VC_snRNAseq_9_Animal4_L-Comb                      73.7  26.3 100.0
# VC_snRNAseq_9_Animal4_LD                          56.8  43.2 100.0
# VC_snRNAseq_9_Animal4_LV1                         53.4  46.6 100.0
# VC_snRNAseq_9_Animal4_LV2                         76.8  23.2 100.0
# VC_snRNAseq-7_AccBasalAP1AP2_-7                   78.6  21.4 100.0
# VC_snRNAseq-7_BasalDorsalAP1_-3A                  89.8  10.2 100.0
# VC_snRNAseq-7_BasalDorsalAP1_-3B                  88.3  11.7 100.0
# VC_snRNAseq-7_BasalVentralAP1_-4                  87.9  12.1 100.0
# VC_snRNAseq-7_BasalVentralAP2_-5                  85.3  14.7 100.0
# VC_snRNAseq-7_LateralDorsalAP1AP2_-8              90.1   9.9 100.0
# VC_snRNAseq-7_LateralVentralAP1_-1                89.3  10.7 100.0
# VC_snRNAseq-7_LateralVentralAP2_-2                91.6   8.4 100.0
# Sum                                               81.0  19.0 100.0


sce$low_lib <- sce$low_lib | qc.lib
sce$low_genes <- sce$low_genes | qc.genes


#### QC plots ####
#### QC plots ####
png(here(plot_dir, "Violin_Subsets_mito.png"), width=1, height=14, units="in", res=300)
plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "high_mito") +
  ggtitle("Mito Percent")  +
    scale_colour_manual(values = c("grey", "red")) +
      coord_flip()
dev.off()

png(here(plot_dir, "Violin_sum_genes.png"), width=10, height=14, units="in", res=300)
plotColData(sce, x = "Sample", y = "sum", colour_by = "low_lib") +
  scale_y_log10() +
  ggtitle("Total UMIs")  +
    scale_colour_manual(values = c("grey", "red")) +
      coord_flip()
dev.off()

png(here(plot_dir, "Violin_detected_genes.png"), width=10, height=14, units="in", res=300)
plotColData(sce, x = "Sample", y = "detected", colour_by = "low_genes") +
  scale_y_log10() +
  ggtitle("Detected genes")  +
    scale_colour_manual(values = c("grey", "red")) +
      coord_flip()
dev.off()




# Mito rate vs n detected features
pdf(height=7.5, width=7.5, here(plot_dir, "QC_scatter_plots.pdf"))
plotColData(sce,
            x = "detected", y = "subsets_Mito_percent",
            colour_by = "discard", point_size = 2.5, point_alpha = 0.5
) + 
    ggtitle("Mito percent") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Detected features vs total count
plotColData(sce,
            x = "sum", y = "detected",
            colour_by = "discard", point_size = 2.5, point_alpha = 0.5
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




table(sce$discard, sce$doubletScore >= 2.75)
#       FALSE  TRUE
# FALSE 20821   447
# TRUE   4107    27


#### Save clean data as HDF5 file  ####
#load(here(processed_dir, "sce_no_empty_droplets.Rdata"))

sce <- sce[, !sce$discard]
dim(sce)
# [1]  21369 107958

save(sce,file=here(processed_dir,"sce_post_qc.rda"))


# sgejobs::job_single('03_droplet_qc', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 03_droplet_qc.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
