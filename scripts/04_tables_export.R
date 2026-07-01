############################################################
# 04_tables_export.R
# Export final manuscript and supplementary tables
############################################################

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(tibble)

options(scipen = 999)

############################################################
# 0. Paths
############################################################

# GitHub-ready path handling:
# - Set PROJECT_DIR to the repository root, or run scripts from the repository root.
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = getwd())
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = FALSE)

DATA_DIR    <- file.path(PROJECT_DIR, "data", "processed")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

############################################################
# Inputs from scripts 01-03
############################################################

IN_TABLE1 <- file.path(RESULTS_DIR, "Table1_baseline_v2.csv")

IN_TABLE2_GLM_MAIN <- file.path(RESULTS_DIR, "Table2_glm_main_v2.csv")
IN_TABLE2_GLM_INT  <- file.path(RESULTS_DIR, "Table2_glm_interaction_v2.csv")

IN_METRICS <- file.path(RESULTS_DIR, "model_metrics_validation_v2.csv")
IN_DELONG  <- file.path(RESULTS_DIR, "delong_auc_comparisons_v2.csv")
IN_AIC     <- file.path(RESULTS_DIR, "model_aic_full_v2.csv")
IN_GLM_LRT <- file.path(RESULTS_DIR, "glm_interaction_lrt_v2.csv")
IN_GAM_ANOVA <- file.path(RESULTS_DIR, "gam_interaction_anova_v2.csv")
IN_GAM_SMOOTHS <- file.path(RESULTS_DIR, "gam_smooth_terms_v2.csv")

IN_FIG_GAM_MAIN <- file.path(RESULTS_DIR, "Figure2_gam_main_v2.png")
IN_FIG_GAM_INT  <- file.path(RESULTS_DIR, "Figure3_gam_interaction_v2.png")

############################################################
# Outputs
############################################################

OUT_TABLE1_FINAL <- file.path(RESULTS_DIR, "Final_Table1_baseline_characteristics_v2.csv")
OUT_TABLE2_FINAL <- file.path(RESULTS_DIR, "Final_Table2_multivariable_glm_main_v2.csv")

OUT_SUPP_TABLE_GLM_INT <- file.path(RESULTS_DIR, "Supplementary_Table_GLM_interaction_v2.csv")
OUT_SUPP_TABLE_PERFORMANCE <- file.path(RESULTS_DIR, "Supplementary_Table_model_performance_v2.csv")
OUT_SUPP_TABLE_DELONG <- file.path(RESULTS_DIR, "Supplementary_Table_DeLong_AUROC_comparisons_v2.csv")
OUT_SUPP_TABLE_INTERACTION_TESTS <- file.path(RESULTS_DIR, "Supplementary_Table_interaction_tests_v2.csv")
OUT_SUPP_TABLE_AIC <- file.path(RESULTS_DIR, "Supplementary_Table_model_AIC_v2.csv")
OUT_SUPP_TABLE_GAM_SMOOTHS <- file.path(RESULTS_DIR, "Supplementary_Table_GAM_smooth_terms_v2.csv")

OUT_MARKDOWN <- file.path(RESULTS_DIR, "final_tables_for_manuscript_v2.md")
OUT_DOCX     <- file.path(RESULTS_DIR, "final_tables_for_manuscript_v2.docx")
OUT_AUDIT    <- file.path(RESULTS_DIR, "tables_export_audit_v2.txt")

############################################################
# 1. Helper functions
############################################################

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop(paste0("Missing required input file: ", path))
  }
}

fmt_num <- function(x, digits = 3) {
  x <- as.numeric(x)
  
  ifelse(
    is.na(x),
    "",
    trimws(format(round(x, digits), nsmall = digits))
  )
}

fmt_auc <- function(auc, low, high) {
  paste0(
    fmt_num(auc, 3),
    " (",
    fmt_num(low, 3),
    "–",
    fmt_num(high, 3),
    ")"
  )
}

fmt_p <- function(p) {
  p <- as.numeric(p)
  
  ifelse(
    is.na(p),
    "",
    ifelse(
      p < 0.001,
      "<0.001",
      trimws(format(round(p, 3), nsmall = 3))
    )
  )
}

fmt_difference <- function(x) {
  x <- as.numeric(x)
  
  ifelse(
    is.na(x),
    "",
    trimws(format(round(x, 4), nsmall = 4))
  )
}

clean_for_export <- function(df) {
  df %>%
    mutate(
      across(
        everything(),
        ~ ifelse(is.na(.x), "", as.character(.x))
      )
    )
}

format_or_table <- function(df) {
  df %>%
    mutate(
      OR = trimws(format(round(as.numeric(OR), 2), nsmall = 2)),
      `95% CI` = as.character(`95% CI`),
      `p-value` = as.character(`p-value`)
    ) %>%
    select(
      Variable,
      OR,
      `95% CI`,
      `p-value`
    ) %>%
    clean_for_export()
}

############################################################
# 2. Check required inputs
############################################################

input_files <- c(
  IN_TABLE1,
  IN_TABLE2_GLM_MAIN,
  IN_TABLE2_GLM_INT,
  IN_METRICS,
  IN_DELONG,
  IN_AIC,
  IN_GLM_LRT,
  IN_GAM_ANOVA,
  IN_GAM_SMOOTHS
)

purrr::walk(input_files, check_file_exists)

############################################################
# 3. Read inputs
############################################################

table1_raw <- read_csv(IN_TABLE1, show_col_types = FALSE)
table2_glm_main_raw <- read_csv(IN_TABLE2_GLM_MAIN, show_col_types = FALSE)
table2_glm_int_raw  <- read_csv(IN_TABLE2_GLM_INT, show_col_types = FALSE)

metrics_raw <- read_csv(IN_METRICS, show_col_types = FALSE)
delong_raw  <- read_csv(IN_DELONG, show_col_types = FALSE)
aic_raw     <- read_csv(IN_AIC, show_col_types = FALSE)
glm_lrt_raw <- read_csv(IN_GLM_LRT, show_col_types = FALSE)
gam_anova_raw <- read_csv(IN_GAM_ANOVA, show_col_types = FALSE)
gam_smooths_raw <- read_csv(IN_GAM_SMOOTHS, show_col_types = FALSE)

############################################################
# 4. Final Table 1
############################################################

table1_final <- table1_raw %>%
  clean_for_export() %>%
  mutate(
    Level = ifelse(Level == "—", "", Level)
  )

write_csv(table1_final, OUT_TABLE1_FINAL)

############################################################
# 5. Final Table 2: primary GLM main-effects model
############################################################

table2_final <- table2_glm_main_raw %>%
  format_or_table()

write_csv(table2_final, OUT_TABLE2_FINAL)

############################################################
# 6. Supplementary table: GLM interaction model
############################################################

supp_glm_interaction <- table2_glm_int_raw %>%
  format_or_table()

write_csv(supp_glm_interaction, OUT_SUPP_TABLE_GLM_INT)

############################################################
# 7. Supplementary model performance table
############################################################

model_label <- function(x) {
  recode(
    x,
    "GLM_main" = "GLM main effects",
    "GLM_interaction" = "GLM with ED LOS × vasopressors interaction",
    "GAM_main" = "GAM main effects",
    "GAM_interaction" = "GAM with ED LOS × vasopressors interaction",
    .default = x
  )
}

supp_performance <- metrics_raw %>%
  mutate(
    Model = model_label(model),
    `Validation N` = as.character(validation_n),
    `Events in validation set` = as.character(validation_events),
    `AUROC (95% CI)` = fmt_auc(auc, auc_ci_low, auc_ci_high),
    `Brier score` = fmt_num(brier, 4),
    `Calibration intercept` = fmt_num(calibration_intercept, 3),
    `Calibration slope` = fmt_num(calibration_slope, 3)
  ) %>%
  select(
    Model,
    `Validation N`,
    `Events in validation set`,
    `AUROC (95% CI)`,
    `Brier score`,
    `Calibration intercept`,
    `Calibration slope`
  ) %>%
  clean_for_export()

write_csv(supp_performance, OUT_SUPP_TABLE_PERFORMANCE)

############################################################
# 8. Supplementary DeLong AUROC comparison table
############################################################

supp_delong <- delong_raw %>%
  mutate(
    Comparison = comparison,
    `Model A` = model_label(model_a),
    `Model B` = model_label(model_b),
    `AUROC model A` = fmt_num(auc_a, 3),
    `AUROC model B` = fmt_num(auc_b, 3),
    `AUROC difference` = fmt_difference(auc_difference_a_minus_b),
    `95% CI for difference` = paste0(
      fmt_difference(ci_low),
      "–",
      fmt_difference(ci_high)
    ),
    `p-value` = fmt_p(p_value)
  ) %>%
  select(
    Comparison,
    `Model A`,
    `Model B`,
    `AUROC model A`,
    `AUROC model B`,
    `AUROC difference`,
    `95% CI for difference`,
    `p-value`
  ) %>%
  clean_for_export()

write_csv(supp_delong, OUT_SUPP_TABLE_DELONG)

############################################################
# 9. Supplementary interaction tests table
############################################################

glm_lrt_clean <- glm_lrt_raw %>%
  filter(!is.na(`Pr(>Chi)`)) %>%
  transmute(
    Test = "GLM likelihood ratio test: main effects vs ED LOS × vasopressors interaction",
    `Test statistic` = fmt_num(Deviance, 3),
    `Degrees of freedom` = fmt_num(Df, 3),
    `p-value` = fmt_p(`Pr(>Chi)`)
  )

gam_anova_clean <- gam_anova_raw %>%
  filter(!is.na(`Pr(>Chi)`)) %>%
  transmute(
    Test = "GAM analysis of deviance: main effects vs ED LOS × vasopressors interaction",
    `Test statistic` = fmt_num(Deviance, 3),
    `Degrees of freedom` = fmt_num(Df, 3),
    `p-value` = fmt_p(`Pr(>Chi)`)
  )

supp_interaction_tests <- bind_rows(
  glm_lrt_clean,
  gam_anova_clean
) %>%
  clean_for_export()

write_csv(supp_interaction_tests, OUT_SUPP_TABLE_INTERACTION_TESTS)

############################################################
# 10. Supplementary AIC table
############################################################

supp_aic <- aic_raw %>%
  mutate(
    Model = model_label(model),
    AIC = fmt_num(AIC, 1)
  ) %>%
  select(
    Model,
    AIC
  ) %>%
  clean_for_export()

write_csv(supp_aic, OUT_SUPP_TABLE_AIC)

############################################################
# 11. Supplementary GAM smooth terms table
############################################################

smooth_p_col <- intersect(
  names(gam_smooths_raw),
  c("p-value", "p.value", "p_value")
)

if (length(smooth_p_col) == 0) {
  gam_smooths_raw <- gam_smooths_raw %>%
    mutate(`p-value` = NA_real_)
  smooth_p_col <- "p-value"
}

supp_gam_smooths <- gam_smooths_raw %>%
  mutate(
    Model = model_label(model),
    `Smooth term` = smooth_term,
    EDF = fmt_num(edf, 3),
    `Reference df` = fmt_num(`Ref.df`, 3),
    `Chi-square` = fmt_num(Chi.sq, 3),
    `p-value` = fmt_p(.data[[smooth_p_col[1]]])
  ) %>%
  select(
    Model,
    `Smooth term`,
    EDF,
    `Reference df`,
    `Chi-square`,
    `p-value`
  ) %>%
  clean_for_export()

write_csv(supp_gam_smooths, OUT_SUPP_TABLE_GAM_SMOOTHS)

############################################################
# 12. Markdown export
############################################################

write_md_table <- function(df, title, con) {
  
  df <- clean_for_export(df)
  
  cat("\n\n", title, "\n", file = con, append = TRUE, sep = "")
  cat(strrep("-", nchar(title)), "\n\n", file = con, append = TRUE, sep = "")
  
  if (nrow(df) == 0) {
    cat("_No rows._\n", file = con, append = TRUE)
    return(invisible(NULL))
  }
  
  header <- paste(names(df), collapse = " | ")
  divider <- paste(rep("---", ncol(df)), collapse = " | ")
  
  cat(header, "\n", file = con, append = TRUE)
  cat(divider, "\n", file = con, append = TRUE)
  
  for (i in seq_len(nrow(df))) {
    row_vals <- unlist(df[i, ], use.names = FALSE)
    row_vals <- as.character(row_vals)
    row_vals <- str_replace_all(row_vals, "\\|", "\\\\|")
    
    cat(
      paste(row_vals, collapse = " | "),
      "\n",
      file = con,
      append = TRUE
    )
  }
}

cat("# Final manuscript tables v2\n", file = OUT_MARKDOWN)

cat(
  "\nGenerated from scripts 01-03 outputs.\n",
  file = OUT_MARKDOWN,
  append = TRUE
)

write_md_table(
  table1_final,
  "Table 1. Baseline characteristics of the study cohort",
  OUT_MARKDOWN
)

cat(
  "\nContinuous variables are reported as mean (SD), except ED LOS, reported as median [IQR]. Categorical variables are reported as n (%). Between-group comparisons used Student's t-test or Mann-Whitney U test for continuous variables, as appropriate, and chi-square tests for categorical variables.\n",
  file = OUT_MARKDOWN,
  append = TRUE
)

write_md_table(
  table2_final,
  "Table 2. Multivariable logistic regression model for ICU mortality",
  OUT_MARKDOWN
)

cat(
  "\nOdds ratios are adjusted for ED LOS, vasopressor use, sex, age, ethnicity, APACHE IV score, and Charlson Comorbidity Index. The primary model is the main-effects GLM.\n",
  file = OUT_MARKDOWN,
  append = TRUE
)

write_md_table(
  supp_performance,
  "Supplementary Table 1. Validation performance of candidate models",
  OUT_MARKDOWN
)

write_md_table(
  supp_delong,
  "Supplementary Table 2. DeLong comparisons of validation AUROC",
  OUT_MARKDOWN
)

write_md_table(
  supp_interaction_tests,
  "Supplementary Table 3. Formal tests for ED LOS × vasopressors interaction",
  OUT_MARKDOWN
)

write_md_table(
  supp_aic,
  "Supplementary Table 4. Full-data model AIC",
  OUT_MARKDOWN
)

write_md_table(
  supp_gam_smooths,
  "Supplementary Table 5. GAM smooth terms",
  OUT_MARKDOWN
)

write_md_table(
  supp_glm_interaction,
  "Supplementary Table 6. GLM with ED LOS × vasopressors interaction",
  OUT_MARKDOWN
)

############################################################
# 13. Optional Word DOCX export
############################################################

docx_available <- requireNamespace("officer", quietly = TRUE) &&
  requireNamespace("flextable", quietly = TRUE)

if (docx_available) {
  
  doc <- officer::read_docx()
  
  add_ft <- function(doc, df, title, footnote = NULL) {
    
    df <- clean_for_export(df)
    
    doc <- officer::body_add_par(doc, title, style = "heading 2")
    
    ft <- flextable::flextable(df)
    ft <- flextable::theme_booktabs(ft)
    ft <- flextable::autofit(ft)
    ft <- flextable::fontsize(ft, size = 9, part = "all")
    ft <- flextable::align(ft, align = "left", part = "all")
    
    doc <- flextable::body_add_flextable(doc, value = ft)
    
    if (!is.null(footnote)) {
      doc <- officer::body_add_par(doc, footnote, style = "Normal")
    }
    
    doc <- officer::body_add_par(doc, "", style = "Normal")
    return(doc)
  }
  
  doc <- officer::body_add_par(
    doc,
    "Final manuscript tables v2",
    style = "heading 1"
  )
  
  doc <- officer::body_add_par(
    doc,
    "Generated from reproducible scripts 01-03.",
    style = "Normal"
  )
  
  doc <- add_ft(
    doc,
    table1_final,
    "Table 1. Baseline characteristics of the study cohort",
    "Continuous variables are reported as mean (SD), except ED LOS, reported as median [IQR]. Categorical variables are reported as n (%). Between-group comparisons used Student's t-test or Mann-Whitney U test for continuous variables, as appropriate, and chi-square tests for categorical variables."
  )
  
  doc <- add_ft(
    doc,
    table2_final,
    "Table 2. Multivariable logistic regression model for ICU mortality",
    "Odds ratios are adjusted for ED LOS, vasopressor use, sex, age, ethnicity, APACHE IV score, and Charlson Comorbidity Index. The primary model is the main-effects GLM."
  )
  
  doc <- add_ft(
    doc,
    supp_performance,
    "Supplementary Table 1. Validation performance of candidate models"
  )
  
  doc <- add_ft(
    doc,
    supp_delong,
    "Supplementary Table 2. DeLong comparisons of validation AUROC"
  )
  
  doc <- add_ft(
    doc,
    supp_interaction_tests,
    "Supplementary Table 3. Formal tests for ED LOS × vasopressors interaction"
  )
  
  doc <- add_ft(
    doc,
    supp_aic,
    "Supplementary Table 4. Full-data model AIC"
  )
  
  doc <- add_ft(
    doc,
    supp_gam_smooths,
    "Supplementary Table 5. GAM smooth terms"
  )
  
  doc <- add_ft(
    doc,
    supp_glm_interaction,
    "Supplementary Table 6. GLM with ED LOS × vasopressors interaction"
  )
  
  if (file.exists(IN_FIG_GAM_MAIN)) {
    doc <- officer::body_add_par(
      doc,
      "Figure 2. GAM main-effect model",
      style = "heading 2"
    )
    doc <- officer::body_add_img(
      doc,
      src = IN_FIG_GAM_MAIN,
      width = 6.5,
      height = 4.8
    )
  }
  
  if (file.exists(IN_FIG_GAM_INT)) {
    doc <- officer::body_add_par(
      doc,
      "Figure 3. GAM interaction model",
      style = "heading 2"
    )
    doc <- officer::body_add_img(
      doc,
      src = IN_FIG_GAM_INT,
      width = 6.5,
      height = 4.8
    )
  }
  

  print(doc, target = OUT_DOCX)
}

############################################################
# 14. Audit report
############################################################

sink(OUT_AUDIT)

cat("\n====================\n")
cat("TABLES EXPORT AUDIT V2\n")
cat("====================\n\n")

cat("====================\n")
cat("INPUT FILES\n")
cat("====================\n")
print(tibble(input_file = input_files, exists = file.exists(input_files)))

cat("\n====================\n")
cat("FINAL TABLE 1\n")
cat("====================\n")
cat("Rows:", nrow(table1_final), "\n")
cat("Columns:", ncol(table1_final), "\n")
print(table1_final, n = nrow(table1_final))

cat("\n====================\n")
cat("FINAL TABLE 2 - GLM MAIN\n")
cat("====================\n")
cat("Rows:", nrow(table2_final), "\n")
cat("Columns:", ncol(table2_final), "\n")
print(table2_final, n = nrow(table2_final))

cat("\n====================\n")
cat("MODEL PERFORMANCE TABLE\n")
cat("====================\n")
print(supp_performance, n = nrow(supp_performance))

cat("\n====================\n")
cat("DELONG TABLE\n")
cat("====================\n")
print(supp_delong, n = nrow(supp_delong))

cat("\n====================\n")
cat("INTERACTION TESTS TABLE\n")
cat("====================\n")
print(supp_interaction_tests, n = nrow(supp_interaction_tests))

cat("\n====================\n")
cat("AIC TABLE\n")
cat("====================\n")
print(supp_aic, n = nrow(supp_aic))

cat("\n====================\n")
cat("GAM SMOOTH TERMS TABLE\n")
cat("====================\n")
print(supp_gam_smooths, n = nrow(supp_gam_smooths))

cat("\n====================\n")
cat("SUPPLEMENTARY GLM INTERACTION TABLE\n")
cat("====================\n")
print(supp_glm_interaction, n = nrow(supp_glm_interaction))

cat("\n====================\n")
cat("DOCX EXPORT\n")
cat("====================\n")
cat("officer/flextable available:", docx_available, "\n")
if (docx_available) {
  cat("DOCX saved to:", OUT_DOCX, "\n")
} else {
  cat("DOCX not created because officer and/or flextable are not installed.\n")
}

cat("\n====================\n")
cat("OUTPUT FILES\n")
cat("====================\n")
cat("Final Table 1:", OUT_TABLE1_FINAL, "\n")
cat("Final Table 2:", OUT_TABLE2_FINAL, "\n")
cat("Supplementary GLM interaction:", OUT_SUPP_TABLE_GLM_INT, "\n")
cat("Supplementary performance:", OUT_SUPP_TABLE_PERFORMANCE, "\n")
cat("Supplementary DeLong:", OUT_SUPP_TABLE_DELONG, "\n")
cat("Supplementary interaction tests:", OUT_SUPP_TABLE_INTERACTION_TESTS, "\n")
cat("Supplementary AIC:", OUT_SUPP_TABLE_AIC, "\n")
cat("Supplementary GAM smooths:", OUT_SUPP_TABLE_GAM_SMOOTHS, "\n")
cat("Markdown:", OUT_MARKDOWN, "\n")
cat("DOCX:", OUT_DOCX, "\n")
cat("Audit:", OUT_AUDIT, "\n")

sink()

############################################################
# 15. Console output
############################################################

cat("\n====================\n")
cat("04 TABLES EXPORT V2 FINISHED\n")
cat("====================\n")
cat("Audit saved to:\n", OUT_AUDIT, "\n")
cat("Markdown saved to:\n", OUT_MARKDOWN, "\n")

if (docx_available) {
  cat("DOCX saved to:\n", OUT_DOCX, "\n")
} else {
  cat("DOCX was not created because officer/flextable are not installed.\n")
}

