library("rsconnect")

#source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "sce_FINAL_human.rda", "initial.R"),
    appName = "BLA_Human",
    account = "libd",
    server = "shinyapps.io"
)
