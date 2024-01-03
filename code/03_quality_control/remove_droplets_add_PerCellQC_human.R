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
load(here("processed-data", "02_build_sce", "sce_human_raw.rda"), verbose = TRUE)

pd <- colData(sce) %>% as.data.frame()

## save directories
plot_dir = here("plots", "03_quality_control","Human","PerCellQC")
processed_dir = here("processed-data","03_quality_control","Human","PerCellQC")

##get droplet score filepaths
droplet_paths <- list.files(here("processed-data","03_quality_control","Human", "droplet_scores"),
                            full.names = TRUE
)

names(droplet_paths) <- gsub("st", "s", gsub("droplet_scores_|.Rdata", "", basename(droplet_paths)))

e.out <- lapply(droplet_paths, function(x) get(load(x)))
e.out

# Read the log file
#log_file_path <- here('code','03_quality_control','logs','emptyDrops_human.txt')
#log_text <- readLines(log_file_path, warn = FALSE)
#log_text

# Find the lines with 'knee_lower' values
#knee_lower_lines <- grep('knee_lower', log_text, value = TRUE)
#knee_lower_lines

# Extract the 'knee_lower' values
#knee_lower_values <- str_extract(knee_lower_lines, "(?<=knee_lower =)\\d+")
#knee_lower_values

#knee_lower <- as.numeric(knee_lower_values)
#names(knee_lower) <- c('BR2327','Br8692','Br9021','Br8337', 'Br5273')
#knee_lower

## This is weird. doesn't work with my logs. Easier to just merge these scripts? - MT
#knee_lower <- map_dbl(logs, ~ parse_number(.x[grepl("knee_lower =", .x)]))
#names(knee_lower) <- gsub("st", "s", map_chr(logs, ~str_sub(.x[grepl("Running Sample: ", .x)], " ", 2)))

knee_lower <- c(220,215,220,248,244)
names(knee_lower) <-c('BR2327','Br8692','Br9021','Br8337', 'Br5273')
knee_lower



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
e.out.all <- do.call("rbind", e.out)[colnames(sce), ]
sce <- sce[, which(e.out.all$FDR <= 0.001)]

dim(sce)
# [1] 36601 25402


#### Compute QC metrics ####
sce <- scuttle::addPerCellQC(
  sce,
  subsets = list(Mito = which(seqnames(sce) == "chrM")),
  # BPPARAM = BiocParallel::MulticoreParam(4)
)




# #save drops removed sce
save(sce,file=here("processed-data","03_quality_control","Human","sce_drops_removed_human.rda"))

load(here("processed-data","03_quality_control","Human","sce_drops_removed_human.rda"))

#### Compute QC metrics ####
location <- rowRanges(sce)
is.mito <- any(seqnames(location)=="chrM")

sce <- scuttle::addPerCellQC(
  sce,
  subsets = list(Mito = is.mito),
  # BPPARAM = BiocParallel::MulticoreParam(4)
)

# # ALTERNATIVELY: using resources in AnnotationHub to retrieve chromosomal
# # locations given the Ensembl IDs; this should yield the same result.
# library(AnnotationHub)
# hs.db <- AnnotationHub()[["AH73881"]]
# chr.loc <- mapIds(hs.db, keys=rownames(sce),
#                   keytype="GENEID", column="SEQNAME")
# is.mito.alt <- which(chr.loc=="chrM")


#### Check for low quality nuc ####
## High mito
# sce$high.mito.sample ## standard name?
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$Sample)

table(sce$high_mito)
# FALSE  TRUE 
# 22060  3342 

table(sce$high_mito, sce$Sample)
#       34ac_scp  35ac_scp  3c-AMYBLA  4c-AMYBLA  5c-AMYBLA
# FALSE     4212     4863      4055      4347      4583
# TRUE       655      649       711       651       676



## low library size
sce$low_lib <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sample)
table(sce$low_lib)
# FALSE  TRUE 
# 24231  1171 

## low detected features
# sce$qc.detected
sce$low_genes <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$Sample)
table(sce$low_genes)
# FALSE  TRUE 
# 23322  2080 



## All low sum are also low detected
table(sce$low_lib, sce$low_genes)
#       FALSE  TRUE
# FALSE 23322   909
# TRUE      0  1171


## Annotate nuc to drop
sce$discard_auto <- sce$high_mito | sce$low_lib | sce$low_genes

table(sce$discard_auto)
# FALSE  TRUE 
# 21268  4134 



# ======= Additionally exclude base on minimum thresholds ========

qc.lib <- sce$sum < 600
qc.genes <- sce$detected < 500

sce$discard_minimum <- qc.lib | qc.genes 

table(sce$discard_minimum)


# combine auto and minimum discards
sce$discard <- sce$discard_auto | sce$discard_minimum

table(sce$discard)
# FALSE  TRUE 
# 21212  4190 

qc_t <- addmargins(table(sce$Sample, sce$discard))

qc_t
#                  FALSE  TRUE   Sum
# Br2327-3c-AMYBLA  4010   756  4766
# Br5273-35ac_scp   4719   793  5512
# Br8331-34ac_scp   3962   905  4867
# Br8692-4c-AMYBLA  4052   946  4998
# Br9021-5c-AMYBLA  4469   790  5259
# Sum              21212  4190 25402



round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)
#                  FALSE  TRUE   Sum
# Br2327-3c-AMYBLA  84.1  15.9 100.0
# Br5273-35ac_scp   85.6  14.4 100.0
# Br8331-34ac_scp   81.4  18.6 100.0
# Br8692-4c-AMYBLA  81.1  18.9 100.0
# Br9021-5c-AMYBLA  85.0  15.0 100.0
# Sum               83.5  16.5 100.0



sce$low_lib <- sce$low_lib | qc.lib
sce$low_genes <- sce$low_genes | qc.genes

#### QC plots ####
pdf(here(plot_dir, "QC_violin_plots.pdf"), width = 10)
## Mito rate
plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "high_mito") +
  ggtitle("Mito Percent") #+
#   facet_wrap(~ sce$round, scales = "free_x", nrow = 1)

# ## low sum
plotColData(sce, x = "Sample", y = "sum", colour_by = "low_lib") +
  scale_y_log10() +
  ggtitle("Total UMIs") #+
#  facet_wrap(~ sce$round, scales = "free_x", nrow = 1)
# +
#   geom_hline(yintercept = 1000) ## hline doesn't work w/ facet_wrap?

# ## low detected
plotColData(sce, x = "Sample", y = "detected", colour_by = "low_genes") +
  scale_y_log10() +
  ggtitle("Detected genes") #+
# geom_hline(yintercept = 500)+
#   facet_wrap(~ sce$round, scales = "free_x", nrow = 1)

dev.off()


pdf(here(plot_dir, "QC_scatter_plots.pdf"), width = 10, height=10)
plotColData(sce,
            x = "detected", y = "subsets_Mito_percent",
            colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5
)

# Detected features vs total count

plotColData(sce,
            x = "sum", y = "detected",
            colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5
)
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
  ggplot(aes(x = reorder(Sample, doubletScore, FUN = median), y = doubletScore)) +
  geom_boxplot() +
  labs(x = "Sample") +
  geom_hline(yintercept = 2.75, color = "red", linetype = "dashed") +
  coord_flip() +
  theme(text = element_text(size = 30)) 

ggsave(dbl_box_plot, filename = here(plot_dir, "doublet_scores_boxplot.png"))

dbl_density_plot <- dbl_df %>%
  ggplot(aes(x = doubletScore)) +
  geom_density() +
  labs(x = "doublet score") +
  facet_grid(Sample ~ .) +
  theme_bw() + 
  theme(text = element_text(size = 30)) 

ggsave(dbl_density_plot, filename = here(plot_dir, "doublet_scores_density.png"), height = 17)

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




## Save
save(sce, file = here::here(processed_dir, "sce_no_empty_droplets.Rdata"))

## Save out sample info for easy access
pd <- colData(sce) %>% as.data.frame()
sample_info <- pd %>%
  group_by(Sample) %>%
  summarize(
    n = length(colnames(sce)),
    n_high_mito = sum(high_mito),
    n_low_sum = sum(low_lib),
    n_low_detect = sum(low_genes),
    n_discard_auto = sum(discard)
  )

write_csv(sample_info, file = here(processed_dir, "sample_info.csv"))

# visualize sample info
#n_boxplot <- sample_info %>%
#  ggplot(aes(x = round, y = n, color = round)) +
#  geom_boxplot(outlier.shape = NA) +
#  geom_point() +
#  ggrepel::geom_text_repel(aes(label = Sample), color = "black") +
  #my_theme +
#  theme(legend.position = "None")

#ggsave(n_boxplot, filename = here(plot_dir, "n_nuclei_boxplot.png"))

#### Save clean data as HDF5 file  ####
load(here(processed_dir, "sce_no_empty_droplets.Rdata"))

sce <- sce[, !sce$discard]
dim(sce)
# [1] 36601 21268

save(sce,file=here("processed-data","03_quality_control","Human","sce_post_qc.rda"))


# sgejobs::job_single('03_droplet_qc', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 03_droplet_qc.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
