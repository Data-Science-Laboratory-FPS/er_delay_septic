############################################################
# run_all.R
# Run complete analysis pipeline
############################################################

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = getwd())
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = FALSE)
Sys.setenv(PROJECT_DIR = PROJECT_DIR)

message("PROJECT_DIR: ", PROJECT_DIR)
message("EICU_DIR: ", Sys.getenv("EICU_DIR", unset = file.path(PROJECT_DIR, "data", "raw", "eicu-crd", "2.0")))

source(file.path(PROJECT_DIR, "scripts", "01_clean_data.R"))
source(file.path(PROJECT_DIR, "scripts", "02_table1.R"))
source(file.path(PROJECT_DIR, "scripts", "03_models.R"))
source(file.path(PROJECT_DIR, "scripts", "04_tables_export.R"))

message("Pipeline finished.")
