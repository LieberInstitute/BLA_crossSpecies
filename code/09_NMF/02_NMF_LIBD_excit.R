library(RcppML)
library(here)

processed_dir <- here("processed-data", "09_NMF")

# get toy brain data
sce.excit.human <- readRDS(here(processed_dir, "sce.excit.human.rds"))
sce.excit.human

# get logcounts
logcounts <- logcounts(sce.excit.human)

# run NMF
print("Starting NMF!")
start_time <- Sys.time()
x<-nmf(logcounts,
       100,
       tol = 1e-06,
       maxit = 1000,
       verbose = TRUE,
       seed = 1512,
       L1 = c(0, 0),
       mask_zeros = FALSE,
       diag = TRUE,
       nonneg = TRUE
)
end_time <- Sys.time()
print(end_time - start_time)
# Time difference of 19.1871 mins

saveRDS(x,file=here(processed_dir, "RcppML_NMF_LIBD_Excit.rds"))
