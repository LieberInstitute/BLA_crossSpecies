library("rsconnect")

#source("../token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "sce_FINAL_macaque.rda", "celltype_colors.rds", "initial.R"),
    appName = "BLA_Macaque",
    account = "libd",
    server = "shinyapps.io"
)
