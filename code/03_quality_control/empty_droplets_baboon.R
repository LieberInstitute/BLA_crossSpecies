library("SingleCellExperiment")
library("DropletUtils")
library("scuttle")
library("tidyverse")
library("here")
library("sessioninfo")


# --------------------- Load sce ---------------------
load(here("processed-data", "02_build_sce", "sce_baboon_raw_ncbi.rda"), verbose = TRUE)
sce


unique_samples <- unique(sce$Sample)

for (sample_run in unique_samples) {
    sce_sample <- sce[, sce$Sample == sample_run]
    
    #--------------------- Run barcodeRanks --------------------- 
    
    bcRanks <- barcodeRanks(sce_sample, fit.bounds = c(10, 1e3))
    
    knee_highest <- metadata(bcRanks)$knee - 200
    message(
        "'First knee point' = ", metadata(bcRanks)$knee, "\n",
        "knee_highest =", knee_highest
    )
    
    knee_higher <- metadata(bcRanks)$knee - 100
    message(
        "'Second knee point' = ", metadata(bcRanks)$knee, "\n",
        "knee_higher =", knee_higher
    )
    
    message(
        "'Third knee point' = ", metadata(bcRanks)$knee, "\n",
        "knee =", metadata(bcRanks)$knee
    )
    
    knee_lower <- metadata(bcRanks)$knee + 100
    message(
        "'Fourth knee point' = ", metadata(bcRanks)$knee, "\n",
        "knee_lower =", knee_lower
    )
    
    knee_lowest <- metadata(bcRanks)$knee + 200
    message(
        "'Fifth knee point' = ", metadata(bcRanks)$knee, "\n",
        "knee_lowest =", knee_lowest
    )
    
    # ---------------------- Run emptyDrops w/ knee + 100 ------------------------
    set.seed(100)
    message("Starting emptyDrops")
    Sys.time()
    e.out <- DropletUtils::emptyDrops(
        sce_sample,
        niters = 25000,
        lower = knee_lower,
        BPPARAM = BiocParallel::MulticoreParam(8)
    )
    message("Done - saving data")
    Sys.time()
    
    if (!dir.exists(here("processed-data", "03_quality_control", "Baboon","droplet_scores"))) dir.create(here("processed-data", "03_quality_control","Baboon", "droplet_scores"))
    save(e.out, file = here("processed-data", "03_quality_control","Baboon", "droplet_scores", paste0("ncbi_droplet_scores_", sample_run, ".Rdata")))
    
    
    # ---------------------------Droplet Plots ----------------------------
    message("QC check")
    FDR_cutoff <- 0.001
    addmargins(table(Signif = e.out$FDR <= FDR_cutoff, Limited = e.out$Limited, useNA = "ifany"))
    
    n_cell_anno <- paste("Non-empty:", sum(e.out$FDR < FDR_cutoff, na.rm = TRUE))
    message(n_cell_anno)
    
    my_theme <- theme_bw() +
        theme(text = element_text(size = 15))
    
    droplet_elbow_plot <- as.data.frame(bcRanks) %>%
        add_column(FDR = e.out$FDR) %>%
        ggplot(aes(x = rank, y = total, color = FDR < FDR_cutoff)) +
        geom_point(alpha = 0.5, size = 1) +
        geom_hline(yintercept = metadata(bcRanks)$knee, linetype = "dotted", color = "gray") +
        annotate("text", x = 10, y = metadata(bcRanks)$knee, label = "Knee", vjust = -1, color = "gray") +
        geom_hline(yintercept = knee_highest, linetype = "dashed") +
        annotate("text", x = 10, y = knee_highest, label = "Knee est 'highest'") +
        geom_hline(yintercept = knee_higher, linetype = "dashed") +
        annotate("text", x = 10, y = knee_higher, label = "Knee est 'higher'") +
        geom_hline(yintercept = knee_lower, linetype = "dashed") +
        annotate("text", x = 10, y = knee_lower, label = "Knee est 'lower'") +
        geom_hline(yintercept = knee_lowest, linetype = "dashed") +
        annotate("text", x = 10, y = knee_lowest, label = "Knee est 'lowest'") +
        scale_x_continuous(trans = "log10") +
        scale_y_continuous(trans = "log10") +
        labs(
            x = "Barcode Rank",
            y = "Total UMIs",
            title = paste("Sample", sample_run),
            subtitle = n_cell_anno,
            color = paste("FDR <", FDR_cutoff)
        ) +
        my_theme +
        theme(legend.position = "bottom")
    
    # droplet_scatter_plot <- as.data.frame(e) %>%
    #   ggplot(aes(x = Total, y = -LogProb, color = FDR < FDR_cutoff)) +
    #   geom_point(alpha = 0.5, size = 1) +
    #   labs(x = "Total UMIs", y = "-Log Probability",
    #        color = paste("FDR <", FDR_cutoff)) +
    #   my_theme+
    #   theme(legend.position = "bottom")
    # # print(droplet_elbow_plot/droplet_scatter_plot)
    # ggsave(droplet_elbow_plot/droplet_scatter_plot, filename = here("plots","03_build_sce", "droplet_qc_png",paste0("droplet_qc_",sample,".png")))
    
    ggsave(droplet_elbow_plot, filename = here("plots", "03_quality_control","Baboon", "droplet_scores", paste0("droplet_qc_", sample_run, ".png")))
}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()