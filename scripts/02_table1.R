############################################################
# 02_table1.R
# Baseline characteristics table using final cleaned cohort v2
############################################################

library(dplyr)
library(readr)
library(tidyr)
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

OUT_TABLE1_CSV  <- file.path(OUT_DIR, "Table1_baseline_v2.csv")
OUT_TABLE1_AUDIT <- file.path(OUT_DIR, "Table1_audit_v2.txt")

############################################################
# 1. Load data
############################################################

df <- read_csv(DATA_PATH, show_col_types = FALSE)

############################################################
# 2. Prepare variables
############################################################

df <- df %>%
  mutate(
    icu_mortality = as.integer(icu_mortality),
    group = case_when(
      icu_mortality == 0 ~ "Survivors",
      icu_mortality == 1 ~ "Non-survivors",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("Female", "Male")),
    vasopressors = factor(
      vasopressors,
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),
    ethnicity = ifelse(
      is.na(ethnicity) | ethnicity == "",
      "Other/Unknown",
      ethnicity
    ),
    ethnicity = str_trim(ethnicity),
    ethnicity = ifelse(
      ethnicity == "",
      "Other/Unknown",
      ethnicity
    ),
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
  filter(!is.na(icu_mortality))

############################################################
# 3. Denominators
############################################################

n_total <- nrow(df)
n_surv  <- sum(df$icu_mortality == 0)
n_nons  <- sum(df$icu_mortality == 1)

col_overall <- paste0("Overall (n = ", format(n_total, big.mark = ","), ")")
col_surv    <- paste0("Survivors (n = ", format(n_surv, big.mark = ","), ")")
col_nons    <- paste0("Non-survivors (n = ", format(n_nons, big.mark = ","), ")")

############################################################
# 4. Formatting helpers
############################################################

fmt_mean_sd <- function(x, digits = 1) {
  paste0(
    format(round(mean(x, na.rm = TRUE), digits), nsmall = digits),
    " (",
    format(round(sd(x, na.rm = TRUE), digits), nsmall = digits),
    ")"
  )
}

fmt_median_iqr <- function(x, digits = 1) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  
  paste0(
    format(round(median(x, na.rm = TRUE), digits), nsmall = digits),
    " [",
    format(round(q1, digits), nsmall = digits),
    "–",
    format(round(q3, digits), nsmall = digits),
    "]"
  )
}

fmt_n_pct <- function(n, denom, digits = 1) {
  pct <- 100 * n / denom
  paste0(
    n,
    " (",
    format(round(pct, digits), nsmall = digits),
    "%)"
  )
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

############################################################
# 5. Statistical test helpers
############################################################

get_ttest_p <- function(var) {
  formula <- as.formula(paste(var, "~ icu_mortality"))
  t.test(formula, data = df)$p.value
}

get_wilcox_p <- function(var) {
  formula <- as.formula(paste(var, "~ icu_mortality"))
  wilcox.test(formula, data = df)$p.value
}

get_cat_p <- function(var) {
  tab <- table(df[[var]], df$icu_mortality)
  
  if (any(tab < 5)) {
    fisher.test(tab)$p.value
  } else {
    chisq.test(tab)$p.value
  }
}

############################################################
# 6. Continuous variables
############################################################
# Main manuscript currently uses mean (SD).
# ED LOS is usually skewed, so we provide median [IQR] for ED LOS.
# If you want all variables as mean (SD), replace fmt_median_iqr with fmt_mean_sd for ED LOS.

continuous_rows <- tibble(
  Variable = c(
    "Age (years)",
    "ED LOS (hours)",
    "APACHE IV score",
    "Charlson Comorbidity Index"
  ),
  Level = "—",
  !!col_overall := c(
    fmt_mean_sd(df$age),
    fmt_median_iqr(df$er_los_hrs),
    fmt_mean_sd(df$apacheiv),
    fmt_mean_sd(df$final_charlson_score)
  ),
  !!col_surv := c(
    fmt_mean_sd(df$age[df$icu_mortality == 0]),
    fmt_median_iqr(df$er_los_hrs[df$icu_mortality == 0]),
    fmt_mean_sd(df$apacheiv[df$icu_mortality == 0]),
    fmt_mean_sd(df$final_charlson_score[df$icu_mortality == 0])
  ),
  !!col_nons := c(
    fmt_mean_sd(df$age[df$icu_mortality == 1]),
    fmt_median_iqr(df$er_los_hrs[df$icu_mortality == 1]),
    fmt_mean_sd(df$apacheiv[df$icu_mortality == 1]),
    fmt_mean_sd(df$final_charlson_score[df$icu_mortality == 1])
  ),
  `p-value` = c(
    fmt_p(get_ttest_p("age")),
    fmt_p(get_wilcox_p("er_los_hrs")),
    fmt_p(get_ttest_p("apacheiv")),
    fmt_p(get_ttest_p("final_charlson_score"))
  )
)

############################################################
# 7. Categorical variable rows
############################################################

make_cat_rows <- function(var, label) {
  
  levels_var <- levels(df[[var]])
  p_value <- fmt_p(get_cat_p(var))
  
  tibble(
    Variable = c(label, rep("", length(levels_var) - 1)),
    Level = levels_var,
    !!col_overall := sapply(
      levels_var,
      function(x) fmt_n_pct(sum(df[[var]] == x, na.rm = TRUE), n_total)
    ),
    !!col_surv := sapply(
      levels_var,
      function(x) fmt_n_pct(
        sum(df[[var]] == x & df$icu_mortality == 0, na.rm = TRUE),
        n_surv
      )
    ),
    !!col_nons := sapply(
      levels_var,
      function(x) fmt_n_pct(
        sum(df[[var]] == x & df$icu_mortality == 1, na.rm = TRUE),
        n_nons
      )
    ),
    `p-value` = c(p_value, rep("", length(levels_var) - 1))
  )
}

sex_rows <- make_cat_rows("sex", "Sex")
eth_rows <- make_cat_rows("ethnicity", "Ethnicity")
vaso_rows <- make_cat_rows("vasopressors", "Vasopressors")

############################################################
# 8. Final Table 1
############################################################

table1 <- bind_rows(
  continuous_rows,
  sex_rows,
  eth_rows,
  vaso_rows
)

write_csv(table1, OUT_TABLE1_CSV)

############################################################
# 9. Audit checks
############################################################

ethnicity_check <- df %>%
  count(ethnicity, icu_mortality, name = "n") %>%
  arrange(ethnicity, icu_mortality)

sex_check <- df %>%
  count(sex, icu_mortality, name = "n") %>%
  arrange(sex, icu_mortality)

vasopressor_check <- df %>%
  count(vasopressors, icu_mortality, name = "n") %>%
  arrange(vasopressors, icu_mortality)

missing_check <- df %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "n_missing"
  )

continuous_summary <- df %>%
  summarise(
    n = n(),
    age_mean = mean(age),
    age_sd = sd(age),
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

sink(OUT_TABLE1_AUDIT)

cat("\n====================\n")
cat("TABLE 1 AUDIT V2\n")
cat("====================\n\n")

cat("Input dataset:\n")
cat(DATA_PATH, "\n\n")

cat("====================\n")
cat("DENOMINATORS\n")
cat("====================\n")
cat("Overall:", n_total, "\n")
cat("Survivors:", n_surv, "\n")
cat("Non-survivors:", n_nons, "\n")

cat("\n====================\n")
cat("MORTALITY RATE\n")
cat("====================\n")
cat("ICU mortality:", round(100 * n_nons / n_total, 2), "%\n")

cat("\n====================\n")
cat("SEX CHECK\n")
cat("====================\n")
print(sex_check)

cat("\n====================\n")
cat("ETHNICITY CHECK\n")
cat("====================\n")
print(ethnicity_check)

cat("\n====================\n")
cat("VASOPRESSOR CHECK\n")
cat("====================\n")
print(vasopressor_check)

cat("\n====================\n")
cat("CONTINUOUS SUMMARY\n")
cat("====================\n")
print(continuous_summary)

cat("\n====================\n")
cat("MISSING CHECK\n")
cat("====================\n")
print(missing_check)

cat("\n====================\n")
cat("TABLE 1\n")
cat("====================\n")
print(table1, n = nrow(table1))

sink()

############################################################
# 10. Console output
############################################################

cat("\n====================\n")
cat("TABLE 1 SCRIPT FINISHED\n")
cat("====================\n")
cat("Table 1 saved to:\n", OUT_TABLE1_CSV, "\n")
cat("Audit saved to:\n", OUT_TABLE1_AUDIT, "\n\n")

