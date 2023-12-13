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
# [1]  21369 134685

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

load(here(processed_dir, "sce_drops_removed.rda"), verbose = TRUE)

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

table(sce$high_mito, sce$Sample)

# SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1
# FALSE                                        4653
# TRUE                                         1836
# 
# SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1
# FALSE                                       4075
# TRUE                                        2213
# 
# SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1
# FALSE                                       5213
# TRUE                                        2731
# 
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1
# FALSE                                       2178
# TRUE                                        4053
# 
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1
# FALSE                                                 6193
# TRUE                                                  1807
# 
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1
# FALSE                                       9353
# TRUE                                         993
# 
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1
# FALSE                                                10905
# TRUE                                                  2159

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



## All low sum are also low detected
table(sce$low_lib, sce$low_genes)
#       FALSE  TRUE
# FALSE 53487  2052
# TRUE      0  2823

## Annotate nuc to drop
sce$discard_auto <- sce$high_mito | sce$low_lib | sce$low_genes

table(sce$discard_auto)
# FALSE  TRUE 
# 40131 18231 



qc_t <- addmargins(table(sce$Sample, sce$discard_auto))
write_csv(data.frame(qc_t), file = here(processed_dir, "discarded_summary.csv"))

qc_t
#                                                      FALSE  TRUE   Sum
# SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1           4559  1930  6489
# SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1            3916  2372  6288
# SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1            5101  2843  7944
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1            2005  4226  6231
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1  5209  2791  8000
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1            9229  1117 10346
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1 10112  2952 13064
# Sum                                                  40131 18231 58362



round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)

# SNL230508VC_AN__baboon_2_BAMY_1_NeuN_10x_L1           70.3  29.7 100.0
# SNL230508VC_AN_baboon_2_BAMY_2_NeuN_10x_L1            62.3  37.7 100.0
# SNL230508VC_AN_baboon_2_LAMY_1_NeuN_10x_L1            64.2  35.8 100.0
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_10x_L1            32.2  67.8 100.0
# SNL230508VC_HA_baboon_5_BAMY_1_NeuN_plus_DAPI_10x_L1  65.1  34.9 100.0
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_10x_L1            89.2  10.8 100.0
# SNL230508VC_HA_baboon_5_LAMY_1_NeuN_plus_DAPI_10x_L1  77.4  22.6 100.0
# Sum                                                   68.8  31.2 100.0


#### QC plots ####

# Mito rate
png(here(plot_dir, "violin_mito_percent_afterSubset.png"), height=7.5, width=7.5, units="in", res=1200)
plotColData(sce, x = "subsets_Mito_percent", y = "Sample", colour_by = "high_mito") +
    ggtitle("Mito Percent") +
    theme(axis.text.x = element_text(angle = 45))
dev.off()

# low sum
png(here(plot_dir, "violin_sum_afterSubset.png"), height=7.5, width=7.5, units="in", res=1200)
plotColData(sce, x = "sum", y = "Sample", colour_by = "low_lib") +
    scale_y_log10() +
    ggtitle("Total UMIs")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# low detected
png(here(plot_dir, "violin_detected_afterSubset.png"), height=7.5, width=7.5, units="in", res=1200)
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
# [1] 21091 40131

save(sce,file=here(processed_dir,"sce_post_qc_baboon.rda"))


# sgejobs::job_single('03_droplet_qc', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 03_droplet_qc.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
