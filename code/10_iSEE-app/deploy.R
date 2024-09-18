library("rsconnect")

#source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "sce_clean.rda", "initial.R"),
    appName = "LS_snRNAseq",
    account = "libd",
    server = "shinyapps.io"
)
