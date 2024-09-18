library("rsconnect")

#source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "sce_FINAL_all_celltypes.rds", "initial.R"),
    appName = "BLA_crossSpecies",
    account = "libd",
    server = "shinyapps.io"
)
