############################################################
# 03_models.R
# GLM/GAM models with fair interaction comparison
# Validation split 70/30 + AUROC/Brier/DeLong
############################################################

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)
library(mgcv)
library(broom)
library(pROC)
library(stringr)

options(scipen = 999)

############################################################
# 0. Paths
############################################################

# GitHub-ready path handling:
# - Set PROJECT_DIR to the repository root, or run scripts from the repository root.
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = getwd())
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = FALSE)

DATA_DIR <- file.path(PROJECT_DIR, "data", "processed")
OUT_DIR  <- file.path(PROJECT_DIR, "results")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DATA_PATH <- file.path(DATA_DIR, "final_df_er_sepsis_clean_v2.csv")

OUT_TABLE2_GLM_MAIN        <- file.path(OUT_DIR, "Table2_glm_main_v2.csv")
OUT_TABLE2_GLM_INTERACTION <- file.path(OUT_DIR, "Table2_glm_interaction_v2.csv")
OUT_METRICS               <- file.path(OUT_DIR, "model_metrics_validation_v2.csv")
OUT_DELONG                <- file.path(OUT_DIR, "delong_auc_comparisons_v2.csv")
OUT_AIC                   <- file.path(OUT_DIR, "model_aic_full_v2.csv")
OUT_GLM_LRT               <- file.path(OUT_DIR, "glm_interaction_lrt_v2.csv")
OUT_GAM_ANOVA             <- file.path(OUT_DIR, "gam_interaction_anova_v2.csv")
OUT_SMOOTH_TERMS          <- file.path(OUT_DIR, "gam_smooth_terms_v2.csv")
OUT_AUDIT                 <- file.path(OUT_DIR, "model_audit_v2.txt")

OUT_FIG_GAM_MAIN          <- file.path(OUT_DIR, "Figure2_gam_main_v2.png")
OUT_FIG_GAM_INTERACTION   <- file.path(OUT_DIR, "Figure3_gam_interaction_v2.png")

OUT_PRED_GAM_MAIN         <- file.path(OUT_DIR, "prediction_grid_gam_main_v2.csv")
OUT_PRED_GAM_INTERACTION  <- file.path(OUT_DIR, "prediction_grid_gam_interaction_v2.csv")

OUT_MODELS_RDS            <- file.path(OUT_DIR, "models_full_v2.rds")

############################################################
# 1. Load data
############################################################

df_raw <- read_csv(DATA_PATH, show_col_types = FALSE)

required_vars <- c(
  "patientunitstayid",
  "icu_mortality",
  "er_los_hrs",
  "vasopressors",
  "age",
  "sex",
  "ethnicity",
  "apacheiv",
  "final_charlson_score"
)

missing_vars <- setdiff(required_vars, names(df_raw))

if (length(missing_vars) > 0) {
  stop(
    paste0(
      "Missing required variables in input dataset: ",
      paste(missing_vars, collapse = ", ")
    )
  )
}

############################################################
# 2. Prepare analysis dataset
############################################################

df <- df_raw %>%
  mutate(
    icu_mortality = as.integer(icu_mortality),
    er_los_hrs = as.numeric(er_los_hrs),
    age = as.numeric(age),
    apacheiv = as.numeric(apacheiv),
    final_charlson_score = as.numeric(final_charlson_score),
    
    vasopressors = case_when(
      as.character(vasopressors) %in% c("1", "Yes", "yes", "TRUE", "True") ~ "Yes",
      as.character(vasopressors) %in% c("0", "No", "no", "FALSE", "False") ~ "No",
      TRUE ~ NA_character_
    ),
    
    sex = as.character(sex),
    ethnicity = as.character(ethnicity),
    ethnicity = ifelse(
      is.na(ethnicity) | str_trim(ethnicity) == "",
      "Other/Unknown",
      ethnicity
    )
  ) %>%
  filter(
    !is.na(patientunitstayid),
    !is.na(icu_mortality),
    icu_mortality %in% c(0, 1),
    !is.na(er_los_hrs),
    !is.na(age),
    !is.na(sex),
    !is.na(ethnicity),
    !is.na(vasopressors),
    !is.na(apacheiv),
    !is.na(final_charlson_score)
  ) %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    vasopressors = factor(vasopressors, levels = c("No", "Yes")),
    ethnicity = factor(
      ethnicity,
      levels = c(
        "Caucasian",
        "African American",
        "Hispanic",
        "Asian",
        "Native American",
        "Other/Unknown"
      )
    )
  ) %>%
  filter(
    !is.na(sex),
    !is.na(vasopressors),
    !is.na(ethnicity)
  )

if (anyDuplicated(df$patientunitstayid) > 0) {
  warning("Duplicated patientunitstayid values detected in analysis dataset.")
}

if (all(df$icu_mortality == 0) | all(df$icu_mortality == 1)) {
  stop("Outcome has only one class; models cannot be fitted.")
}

############################################################
# 3. Stratified 70/30 split
############################################################

set.seed(20260619)

train_ids <- df %>%
  group_by(icu_mortality) %>%
  slice_sample(prop = 0.70) %>%
  ungroup() %>%
  pull(patientunitstayid)

train_df <- df %>%
  filter(patientunitstayid %in% train_ids)

test_df <- df %>%
  filter(!patientunitstayid %in% train_ids)

split_counts <- bind_rows(
  train_df %>%
    count(icu_mortality, name = "n") %>%
    mutate(dataset = "Training"),
  test_df %>%
    count(icu_mortality, name = "n") %>%
    mutate(dataset = "Validation")
) %>%
  relocate(dataset)

############################################################
# 4. Model formulas
############################################################

formula_glm_main <- icu_mortality ~
  er_los_hrs +
  vasopressors +
  sex +
  age +
  ethnicity +
  apacheiv +
  final_charlson_score

formula_glm_interaction <- icu_mortality ~
  er_los_hrs * vasopressors +
  sex +
  age +
  ethnicity +
  apacheiv +
  final_charlson_score

formula_gam_main <- icu_mortality ~
  s(er_los_hrs, k = 5) +
  vasopressors +
  sex +
  age +
  ethnicity +
  apacheiv +
  final_charlson_score

formula_gam_interaction <- icu_mortality ~
  s(er_los_hrs, by = vasopressors, k = 5) +
  vasopressors +
  sex +
  age +
  ethnicity +
  apacheiv +
  final_charlson_score

############################################################
# 5. Fit full-data models for inference, tables and figures
############################################################

glm_main_full <- glm(
  formula_glm_main,
  data = df,
  family = binomial()
)

glm_interaction_full <- glm(
  formula_glm_interaction,
  data = df,
  family = binomial()
)

gam_main_full <- bam(
  formula_gam_main,
  data = df,
  family = binomial(),
  method = "fREML",
  select = TRUE
)

gam_interaction_full <- bam(
  formula_gam_interaction,
  data = df,
  family = binomial(),
  method = "fREML",
  select = TRUE
)

############################################################
# 6. Fit training models for validation performance
############################################################

glm_main_train <- glm(
  formula_glm_main,
  data = train_df,
  family = binomial()
)

glm_interaction_train <- glm(
  formula_glm_interaction,
  data = train_df,
  family = binomial()
)

gam_main_train <- bam(
  formula_gam_main,
  data = train_df,
  family = binomial(),
  method = "fREML",
  select = TRUE
)

gam_interaction_train <- bam(
  formula_gam_interaction,
  data = train_df,
  family = binomial(),
  method = "fREML",
  select = TRUE
)

############################################################
# 7. Helper functions
############################################################

fmt_num <- function(x, digits = 2) {
  format(round(x, digits), nsmall = digits)
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "",
    ifelse(
      p < 0.001,
      "<0.001",
      format(round(p, 3), nsmall = 3)
    )
  )
}

brier_score <- function(obs, pred) {
  mean((pred - obs)^2)
}

clamp_prob <- function(p, eps = 1e-6) {
  pmin(pmax(p, eps), 1 - eps)
}

calibration_stats <- function(obs, pred) {
  pred_c <- clamp_prob(pred)
  logit_pred <- qlogis(pred_c)
  
  cal_model <- glm(obs ~ logit_pred, family = binomial())
  
  tibble(
    calibration_intercept = unname(coef(cal_model)[1]),
    calibration_slope = unname(coef(cal_model)[2])
  )
}

make_roc <- function(obs, pred) {
  roc(
    response = obs,
    predictor = pred,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
}

model_metric_row <- function(model_name, obs, pred) {
  roc_obj <- make_roc(obs, pred)
  auc_ci <- ci.auc(roc_obj)
  cal <- calibration_stats(obs, pred)
  
  tibble(
    model = model_name,
    validation_n = length(obs),
    validation_events = sum(obs == 1),
    auc = as.numeric(auc(roc_obj)),
    auc_ci_low = as.numeric(auc_ci[1]),
    auc_ci_high = as.numeric(auc_ci[3]),
    brier = brier_score(obs, pred),
    calibration_intercept = cal$calibration_intercept,
    calibration_slope = cal$calibration_slope
  )
}

safe_predict_response <- function(model, newdata) {
  as.numeric(predict(model, newdata = newdata, type = "response"))
}

############################################################
# 8. Validation predictions
############################################################

obs_test <- test_df$icu_mortality

pred_glm_main <- safe_predict_response(glm_main_train, test_df)
pred_glm_interaction <- safe_predict_response(glm_interaction_train, test_df)
pred_gam_main <- safe_predict_response(gam_main_train, test_df)
pred_gam_interaction <- safe_predict_response(gam_interaction_train, test_df)

valid_idx <- complete.cases(
  obs_test,
  pred_glm_main,
  pred_glm_interaction,
  pred_gam_main,
  pred_gam_interaction
)

obs_test <- obs_test[valid_idx]
pred_glm_main <- pred_glm_main[valid_idx]
pred_glm_interaction <- pred_glm_interaction[valid_idx]
pred_gam_main <- pred_gam_main[valid_idx]
pred_gam_interaction <- pred_gam_interaction[valid_idx]

roc_glm_main <- make_roc(obs_test, pred_glm_main)
roc_glm_interaction <- make_roc(obs_test, pred_glm_interaction)
roc_gam_main <- make_roc(obs_test, pred_gam_main)
roc_gam_interaction <- make_roc(obs_test, pred_gam_interaction)

metrics <- bind_rows(
  model_metric_row("GLM_main", obs_test, pred_glm_main),
  model_metric_row("GLM_interaction", obs_test, pred_glm_interaction),
  model_metric_row("GAM_main", obs_test, pred_gam_main),
  model_metric_row("GAM_interaction", obs_test, pred_gam_interaction)
)

write_csv(metrics, OUT_METRICS)

############################################################
# 9. DeLong AUROC comparisons
############################################################

delong_compare <- function(name_a, roc_a, name_b, roc_b) {
  test <- roc.test(
    roc_a,
    roc_b,
    method = "delong",
    paired = TRUE
  )
  
  ci_vals <- if (!is.null(test$conf.int)) {
    as.numeric(test$conf.int)
  } else {
    c(NA_real_, NA_real_, NA_real_)
  }
  
  tibble(
    comparison = paste(name_a, "vs", name_b),
    model_a = name_a,
    model_b = name_b,
    auc_a = as.numeric(auc(roc_a)),
    auc_b = as.numeric(auc(roc_b)),
    auc_difference_a_minus_b = as.numeric(auc(roc_a)) - as.numeric(auc(roc_b)),
    ci_low = ci_vals[1],
    ci_high = ci_vals[length(ci_vals)],
    p_value = as.numeric(test$p.value)
  )
}

delong_results <- bind_rows(
  delong_compare(
    "GLM_main",
    roc_glm_main,
    "GLM_interaction",
    roc_glm_interaction
  ),
  delong_compare(
    "GAM_main",
    roc_gam_main,
    "GAM_interaction",
    roc_gam_interaction
  ),
  delong_compare(
    "GLM_main",
    roc_glm_main,
    "GAM_main",
    roc_gam_main
  ),
  delong_compare(
    "GLM_interaction",
    roc_glm_interaction,
    "GAM_interaction",
    roc_gam_interaction
  )
) %>%
  mutate(
    p_value_formatted = fmt_p(p_value)
  )

write_csv(delong_results, OUT_DELONG)

############################################################
# 10. Full-data model comparisons
############################################################

model_aic <- tibble(
  model = c(
    "GLM_main",
    "GLM_interaction",
    "GAM_main",
    "GAM_interaction"
  ),
  AIC = c(
    AIC(glm_main_full),
    AIC(glm_interaction_full),
    AIC(gam_main_full),
    AIC(gam_interaction_full)
  )
) %>%
  arrange(AIC)

write_csv(model_aic, OUT_AIC)

glm_lrt <- anova(
  glm_main_full,
  glm_interaction_full,
  test = "Chisq"
) %>%
  as.data.frame() %>%
  rownames_to_column("model_step")

write_csv(glm_lrt, OUT_GLM_LRT)

gam_anova <- anova(
  gam_main_full,
  gam_interaction_full,
  test = "Chisq"
) %>%
  as.data.frame() %>%
  rownames_to_column("model_step")

write_csv(gam_anova, OUT_GAM_ANOVA)

############################################################
# 11. Table 2: GLM odds ratios
############################################################

label_terms <- function(term) {
  recode(
    term,
    "er_los_hrs" = "ED LOS (per hour)",
    "vasopressorsYes" = "Vasopressors (Yes vs No)",
    "sexMale" = "Male sex",
    "age" = "Age (per year)",
    "apacheiv" = "APACHE IV score (per point)",
    "final_charlson_score" = "Charlson Comorbidity Index (per point)",
    "ethnicityAfrican American" = "Ethnicity: African American vs Caucasian",
    "ethnicityHispanic" = "Ethnicity: Hispanic vs Caucasian",
    "ethnicityAsian" = "Ethnicity: Asian vs Caucasian",
    "ethnicityNative American" = "Ethnicity: Native American vs Caucasian",
    "ethnicityOther/Unknown" = "Ethnicity: Other/Unknown vs Caucasian",
    "er_los_hrs:vasopressorsYes" = "ED LOS × vasopressors",
    .default = term
  )
}

format_glm_table <- function(model, model_name) {
  tidy(
    model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      model = model_name,
      Variable = label_terms(term),
      OR = fmt_num(estimate, 2),
      `95% CI` = paste0(
        fmt_num(conf.low, 2),
        "–",
        fmt_num(conf.high, 2)
      ),
      `p-value` = fmt_p(p.value)
    ) %>%
    select(
      model,
      term,
      Variable,
      OR,
      `95% CI`,
      `p-value`,
      estimate,
      conf.low,
      conf.high,
      p.value
    )
}

table2_glm_main <- format_glm_table(
  glm_main_full,
  "GLM_main"
)

table2_glm_interaction <- format_glm_table(
  glm_interaction_full,
  "GLM_interaction"
)

write_csv(table2_glm_main, OUT_TABLE2_GLM_MAIN)
write_csv(table2_glm_interaction, OUT_TABLE2_GLM_INTERACTION)

############################################################
# 12. GAM smooth term summaries
############################################################

extract_smooth_table <- function(model, model_name) {
  s_tab <- summary(model)$s.table
  
  if (is.null(s_tab)) {
    return(tibble())
  }
  
  as.data.frame(s_tab) %>%
    rownames_to_column("smooth_term") %>%
    as_tibble() %>%
    mutate(model = model_name) %>%
    relocate(model, smooth_term)
}

smooth_terms <- bind_rows(
  extract_smooth_table(gam_main_full, "GAM_main"),
  extract_smooth_table(gam_interaction_full, "GAM_interaction")
)

write_csv(smooth_terms, OUT_SMOOTH_TERMS)

############################################################
# 13. Prediction grids for figures
############################################################

reference_values <- tibble(
  sex = factor("Female", levels = levels(df$sex)),
  ethnicity = factor("Caucasian", levels = levels(df$ethnicity)),
  age = mean(df$age, na.rm = TRUE),
  apacheiv = mean(df$apacheiv, na.rm = TRUE),
  final_charlson_score = mean(df$final_charlson_score, na.rm = TRUE)
)

ed_seq <- seq(
  min(df$er_los_hrs, na.rm = TRUE),
  max(df$er_los_hrs, na.rm = TRUE),
  length.out = 200
)

############################################################
# Figure 2: GAM main effect
############################################################

grid_gam_main <- expand.grid(
  er_los_hrs = ed_seq,
  vasopressors = factor("No", levels = levels(df$vasopressors))
) %>%
  as_tibble() %>%
  mutate(
    sex = factor("Female", levels = levels(df$sex)),
    ethnicity = factor("Caucasian", levels = levels(df$ethnicity)),
    age = mean(df$age, na.rm = TRUE),
    apacheiv = mean(df$apacheiv, na.rm = TRUE),
    final_charlson_score = mean(df$final_charlson_score, na.rm = TRUE)
  )

pred_main_link <- predict(
  gam_main_full,
  newdata = grid_gam_main,
  type = "link",
  se.fit = TRUE
)

grid_gam_main <- grid_gam_main %>%
  mutate(
    fit = plogis(as.numeric(pred_main_link$fit)),
    lower = plogis(
      as.numeric(pred_main_link$fit) -
        1.96 * as.numeric(pred_main_link$se.fit)
    ),
    upper = plogis(
      as.numeric(pred_main_link$fit) +
        1.96 * as.numeric(pred_main_link$se.fit)
    )
  )

write_csv(grid_gam_main, OUT_PRED_GAM_MAIN)

fig_gam_main <- ggplot(
  grid_gam_main,
  aes(x = er_los_hrs, y = fit)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.20
  ) +
  geom_line(linewidth = 1.1) +
  labs(
    x = "ED length of stay (hours)",
    y = "Predicted probability of ICU mortality"
  ) +
  theme_bw()

ggsave(
  OUT_FIG_GAM_MAIN,
  fig_gam_main,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

############################################################
# Figure 3: GAM interaction by vasopressors
############################################################

grid_gam_interaction <- expand.grid(
  er_los_hrs = ed_seq,
  vasopressors = factor(levels(df$vasopressors), levels = levels(df$vasopressors))
) %>%
  as_tibble() %>%
  mutate(
    sex = factor("Female", levels = levels(df$sex)),
    ethnicity = factor("Caucasian", levels = levels(df$ethnicity)),
    age = mean(df$age, na.rm = TRUE),
    apacheiv = mean(df$apacheiv, na.rm = TRUE),
    final_charlson_score = mean(df$final_charlson_score, na.rm = TRUE)
  )

pred_int_link <- predict(
  gam_interaction_full,
  newdata = grid_gam_interaction,
  type = "link",
  se.fit = TRUE
)

grid_gam_interaction <- grid_gam_interaction %>%
  mutate(
    fit = plogis(as.numeric(pred_int_link$fit)),
    lower = plogis(
      as.numeric(pred_int_link$fit) -
        1.96 * as.numeric(pred_int_link$se.fit)
    ),
    upper = plogis(
      as.numeric(pred_int_link$fit) +
        1.96 * as.numeric(pred_int_link$se.fit)
    )
  )

write_csv(grid_gam_interaction, OUT_PRED_GAM_INTERACTION)

grid_gam_interaction <- grid_gam_interaction %>%
  mutate(
    vasopressors_label = case_when(
      vasopressors == "No" ~ "No vasopressors",
      vasopressors == "Yes" ~ "Vasopressors",
      TRUE ~ as.character(vasopressors)
    ),
    vasopressors_label = factor(
      vasopressors_label,
      levels = c("No vasopressors", "Vasopressors")
    )
  )

fig_gam_interaction <- ggplot(
  grid_gam_interaction,
  aes(
    x = er_los_hrs,
    y = fit,
    color = vasopressors_label,
    fill = vasopressors_label,
    linetype = vasopressors_label
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 1.1) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(
    x = "ED length of stay (hours)",
    y = "Predicted probability of ICU mortality",
    color = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  theme_bw()
ggsave(
  OUT_FIG_GAM_INTERACTION,
  fig_gam_interaction,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

############################################################
# 14. Save full-data models
############################################################

saveRDS(
  list(
    glm_main_full = glm_main_full,
    glm_interaction_full = glm_interaction_full,
    gam_main_full = gam_main_full,
    gam_interaction_full = gam_interaction_full,
    glm_main_train = glm_main_train,
    glm_interaction_train = glm_interaction_train,
    gam_main_train = gam_main_train,
    gam_interaction_train = gam_interaction_train
  ),
  OUT_MODELS_RDS
)

############################################################
# 15. Audit report
############################################################

sink(OUT_AUDIT)

cat("\n====================\n")
cat("MODEL AUDIT V2\n")
cat("====================\n\n")

cat("Input dataset:\n")
cat(DATA_PATH, "\n\n")

cat("====================\n")
cat("ANALYSIS DATASET\n")
cat("====================\n")
cat("N:", nrow(df), "\n")
cat("Events:", sum(df$icu_mortality == 1), "\n")
cat("Non-events:", sum(df$icu_mortality == 0), "\n")
cat("ICU mortality:", round(100 * mean(df$icu_mortality), 2), "%\n")

cat("\n====================\n")
cat("VARIABLE SUMMARIES\n")
cat("====================\n")
print(
  df %>%
    summarise(
      ed_los_mean = mean(er_los_hrs),
      ed_los_sd = sd(er_los_hrs),
      ed_los_median = median(er_los_hrs),
      ed_los_q1 = quantile(er_los_hrs, 0.25),
      ed_los_q3 = quantile(er_los_hrs, 0.75),
      apache_mean = mean(apacheiv),
      apache_sd = sd(apacheiv),
      charlson_mean = mean(final_charlson_score),
      charlson_sd = sd(final_charlson_score)
    )
)

cat("\n====================\n")
cat("SPLIT COUNTS\n")
cat("====================\n")
print(split_counts)

cat("\n====================\n")
cat("MODEL FORMULAS\n")
cat("====================\n")
cat("\nGLM main:\n")
print(formula_glm_main)
cat("\nGLM interaction:\n")
print(formula_glm_interaction)
cat("\nGAM main:\n")
print(formula_gam_main)
cat("\nGAM interaction:\n")
print(formula_gam_interaction)

cat("\n====================\n")
cat("VALIDATION METRICS\n")
cat("====================\n")
print(metrics)

cat("\n====================\n")
cat("DELONG AUROC COMPARISONS\n")
cat("====================\n")
print(delong_results)

cat("\n====================\n")
cat("FULL-DATA MODEL AIC\n")
cat("====================\n")
print(model_aic)

cat("\n====================\n")
cat("GLM INTERACTION LIKELIHOOD RATIO TEST\n")
cat("====================\n")
print(glm_lrt)

cat("\n====================\n")
cat("GAM INTERACTION ANOVA\n")
cat("====================\n")
print(gam_anova)

cat("\n====================\n")
cat("TABLE 2 - GLM MAIN\n")
cat("====================\n")
print(table2_glm_main, n = nrow(table2_glm_main))

cat("\n====================\n")
cat("TABLE 2 - GLM INTERACTION\n")
cat("====================\n")
print(table2_glm_interaction, n = nrow(table2_glm_interaction))

cat("\n====================\n")
cat("GAM SMOOTH TERMS\n")
cat("====================\n")
print(smooth_terms, n = nrow(smooth_terms))

cat("\n====================\n")
cat("OUTPUT FILES\n")
cat("====================\n")
cat("Table 2 GLM main:", OUT_TABLE2_GLM_MAIN, "\n")
cat("Table 2 GLM interaction:", OUT_TABLE2_GLM_INTERACTION, "\n")
cat("Validation metrics:", OUT_METRICS, "\n")
cat("DeLong comparisons:", OUT_DELONG, "\n")
cat("AIC:", OUT_AIC, "\n")
cat("GLM LRT:", OUT_GLM_LRT, "\n")
cat("GAM ANOVA:", OUT_GAM_ANOVA, "\n")
cat("Smooth terms:", OUT_SMOOTH_TERMS, "\n")
cat("Figure GAM main:", OUT_FIG_GAM_MAIN, "\n")
cat("Figure GAM interaction:", OUT_FIG_GAM_INTERACTION, "\n")
cat("Models RDS:", OUT_MODELS_RDS, "\n")

sink()

############################################################
# 16. Console output
############################################################

cat("\n====================\n")
cat("03 MODELS V2 FINISHED\n")
cat("====================\n")
cat("Audit saved to:\n", OUT_AUDIT, "\n")
cat("Validation metrics saved to:\n", OUT_METRICS, "\n")
cat("DeLong comparisons saved to:\n", OUT_DELONG, "\n")
cat("Table 2 GLM main saved to:\n", OUT_TABLE2_GLM_MAIN, "\n")
cat("Table 2 GLM interaction saved to:\n", OUT_TABLE2_GLM_INTERACTION, "\n\n")

