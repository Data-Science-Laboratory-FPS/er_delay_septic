############################################################
# 01_clean_data.R
# Build analytic cohort and audit checks
############################################################

library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(tidyr)

options(scipen = 999)

############################################################
# 0. Paths
############################################################

# GitHub-ready path handling:
# - Set PROJECT_DIR to the repository root, or run scripts from the repository root.
# - Set EICU_DIR to the local eICU-CRD 2.0 folder.
#   Example:
#   Sys.setenv(EICU_DIR = "/path/to/eicu-crd/2.0")
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = getwd())
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = FALSE)

EICU_DIR <- Sys.getenv(
  "EICU_DIR",
  unset = file.path(PROJECT_DIR, "data", "raw", "eicu-crd", "2.0")
)
EICU_DIR <- normalizePath(EICU_DIR, mustWork = FALSE)

OUT_DIR     <- file.path(PROJECT_DIR, "data", "processed")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_CSV_COHORT <- file.path(OUT_DIR, "final_df_er_sepsis_clean_v2.csv")
OUT_FLOWCHART  <- file.path(OUT_DIR, "flowchart_numbers_v2.csv")
OUT_DX_COUNTS  <- file.path(OUT_DIR, "diagnosis_counts_v2.csv")
OUT_AUDIT_TXT  <- file.path(OUT_DIR, "audit_summary_v2.txt")
OUT_MISSING    <- file.path(OUT_DIR, "missing_summary_v2.csv")
OUT_CHARLSON_AUDIT <- file.path(RESULTS_DIR, "charlson_audit_v2.txt")
OUT_PASTHISTORY_PATHS <- file.path(RESULTS_DIR, "pastHistory_unique_paths_v2.csv")

############################################################
# 1. Options
############################################################

# We use the non-age-adjusted Charlson score because age is included
# as an independent covariate in the multivariable models.
USE_AGE_ADJUSTED_CHARLSON <- FALSE

# Keep original approach for now: unitstaytype == "admit".
# We audit unitvisitnumber but do not force unitvisitnumber == 1 by default.
FILTER_UNITVISITNUMBER_1 <- FALSE

# Export a dedicated Charlson/pastHistory audit from the same script.
# This replaces the previous standalone 01b_audit_charlson.R script.
RUN_CHARLSON_AUDIT <- TRUE

############################################################
# 2. Load raw tables
############################################################

patient_raw <- read_csv(file.path(EICU_DIR, "patient.csv.gz"), show_col_types = FALSE)
treatment   <- read_csv(file.path(EICU_DIR, "treatment.csv.gz"), show_col_types = FALSE)
apache_raw  <- read_csv(file.path(EICU_DIR, "apachePatientResult.csv.gz"), show_col_types = FALSE)
pasthistory <- read_csv(file.path(EICU_DIR, "pastHistory.csv.gz"), show_col_types = FALSE)

############################################################
# 3. Helper for cohort counts
############################################################

flow <- tibble(step = character(), N = integer())

add_flow <- function(step_name, data) {
  flow <<- bind_rows(flow, tibble(step = step_name, N = nrow(data)))
}

############################################################
# 4. Initial cohort
############################################################

p0 <- patient_raw
add_flow("Initial ICU stays", p0)

############################################################
# 5. First ICU admission / admitted unit stay
############################################################

p1 <- p0 %>%
  filter(unitstaytype == "admit")

if (FILTER_UNITVISITNUMBER_1 && "unitvisitnumber" %in% names(p1)) {
  p1 <- p1 %>%
    filter(unitvisitnumber == 1)
}

add_flow("unitstaytype == admit", p1)

############################################################
# 6. Sepsis / severe infection diagnosis
############################################################

infection_dx <- c(
  "Cellulitis and localized soft tissue infections, surgery for",
  "Abscess/infection-cranial, surgery for",
  "Renal infection/abscess",
  "Cellulitis and localized soft tissue infections",
  "Thoracotomy for thoracic/respiratory infection",
  "Infection/abscess, other surgery for",
  "Abscess, neurologic",
  "ARDS-adult respiratory distress syndrome, non-cardiogenic pulmonary edema",
  "Arthritis, septic",
  "Cholecystectomy/cholangitis, surgery for (gallbladder removal)",
  "Complications of previous GI surgery; surgery for (anastomotic leak, bleeding, abscess, infection, dehiscence, etc.)",
  "Complications of previous open heart surgery (i.e. bleeding, infection etc.)",
  "Complications of previous open-heart surgery, surgery for (i.e. bleeding, infection, mediastinal rewiring,leaking aortic graft etc.)",
  "Diverticular disease, surgery for",
  "Endocarditis",
  "Fistula/abscess, surgery for (not inflammatory bowel disease)",
  "GI abscess/cyst",
  "GI Abscess/cyst-primary, surgery for",
  "GI obstruction, surgery for (including lysis of adhesions)",
  "GI perforation/rupture",
  "GI perforation/rupture, surgery for",
  "Inflammatory bowel disease, surgery for",
  "Meningitis",
  "Myositis, viral",
  "Pancreatitis, surgery for",
  "Peritonitis",
  "Peritonitis, surgery for",
  "Pneumonia, aspiration",
  "Pneumonia, bacterial",
  "Pneumonia, fungal",
  "Pneumonia, other",
  "Pneumonia, parasitic (i.e., Pneumocystic pneumonia)",
  "Pneumonia, viral"
)

p2 <- p1 %>%
  filter(
    str_detect(apacheadmissiondx, regex("sepsis", ignore_case = TRUE)) |
      apacheadmissiondx %in% infection_dx
  )

add_flow("Sepsis or severe infection diagnosis", p2)

############################################################
# 7. ED admission
############################################################

p3 <- p2 %>%
  filter(hospitaladmitsource == "Emergency Department")

add_flow("Admitted from ED", p3)

############################################################
# 8. ED LOS calculation
# Important: do not use abs().
# In eICU, hospitaladmitoffset is usually negative when hospital
# admission occurred before ICU admission.
############################################################

offset_audit_before <- p3 %>%
  summarise(
    n = n(),
    missing_offset = sum(is.na(hospitaladmitoffset)),
    offset_positive = sum(as.numeric(hospitaladmitoffset) > 0, na.rm = TRUE),
    offset_zero = sum(as.numeric(hospitaladmitoffset) == 0, na.rm = TRUE),
    offset_negative = sum(as.numeric(hospitaladmitoffset) < 0, na.rm = TRUE),
    min_offset = min(as.numeric(hospitaladmitoffset), na.rm = TRUE),
    median_offset = median(as.numeric(hospitaladmitoffset), na.rm = TRUE),
    max_offset = max(as.numeric(hospitaladmitoffset), na.rm = TRUE)
  )

p4 <- p3 %>%
  mutate(
    hospitaladmitoffset = as.numeric(hospitaladmitoffset),
    er_los_hrs = -hospitaladmitoffset / 60
  ) %>%
  filter(
    !is.na(er_los_hrs),
    er_los_hrs >= 0,
    er_los_hrs <= 24
  )

add_flow("ED LOS between 0 and 24h", p4)

############################################################
# 9. Known sex
############################################################

p5 <- p4 %>%
  rename(sex = gender) %>%
  filter(sex %in% c("Male", "Female"))

add_flow("Known sex", p5)

############################################################
# 10. Age >= 16
############################################################

p6 <- p5 %>%
  mutate(
    age = ifelse(age == "> 89", "90", age),
    age = suppressWarnings(as.numeric(age))
  ) %>%
  filter(!is.na(age), age >= 16)

add_flow("Age >= 16", p6)

############################################################
# 11. Diagnosis harmonization
############################################################

p6 <- p6 %>%
  mutate(
    apacheadmissiondx = ifelse(
      apacheadmissiondx %in% c("Pneumonia, bacterial", "Sepsis, pulmonary"),
      "Sepsis, pneumonia",
      apacheadmissiondx
    )
  )

############################################################
# 12. Outcome
############################################################

patient_cohort <- p6 %>%
  mutate(
    icu_mortality = case_when(
      unitdischargestatus == "Expired" ~ 1L,
      unitdischargestatus == "Alive" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

############################################################
# 13. Vasopressors in first 24h after ICU admission
############################################################

vasopressor_terms <- c(
  "norepinephrine",
  "epinephrine",
  "dopamine",
  "phenylephrine",
  "vasopressin"
)

vaso_24h <- treatment %>%
  filter(
    treatmentoffset >= 0,
    treatmentoffset <= 1440,
    str_detect(
      tolower(treatmentstring),
      paste(vasopressor_terms, collapse = "|")
    )
  ) %>%
  distinct(patientunitstayid) %>%
  mutate(vasopressors = 1L)

patient_cohort <- patient_cohort %>%
  left_join(vaso_24h, by = "patientunitstayid") %>%
  mutate(vasopressors = coalesce(vasopressors, 0L))

############################################################
# 14. APACHE IV score
############################################################

apache_dup_audit <- apache_raw %>%
  count(patientunitstayid, name = "n_rows") %>%
  count(n_rows, name = "n_patients")

if ("apacheversion" %in% names(apache_raw)) {
  apache_version_audit <- apache_raw %>%
    count(apacheversion, sort = TRUE)
} else {
  apache_version_audit <- tibble(note = "No apacheversion column found")
}

apache_for_join <- apache_raw

if ("apacheversion" %in% names(apache_for_join)) {
  apache_for_join <- apache_for_join %>%
    filter(is.na(apacheversion) | apacheversion == "IV")
}

apache_iv <- apache_for_join %>%
  select(patientunitstayid, apachescore) %>%
  arrange(patientunitstayid) %>%
  distinct(patientunitstayid, .keep_all = TRUE) %>%
  rename(apacheiv = apachescore)

patient_cohort <- patient_cohort %>%
  left_join(apache_iv, by = "patientunitstayid")

############################################################
# 15. Charlson Comorbidity Index - robust version
############################################################

cci_patterns <- list(
  mets6 = "Metastases/",
  aids6 = "AIDS/AIDS$|HIV \\(only\\)/HIV positive",
  liver3 = "Cirrhosis/(UGI bleeding|varices|coma|jaundice|ascites|encephalopathy)",
  stroke2 = "Neurologic/Strokes/",
  renal2 = "Renal\\s+\\(R\\)/Renal (Insufficiency|Failure)",
  dm1 = "Endocrine \\(R\\)/(Insulin Dependent Diabetes|Non-Insulin Dependent Diabetes)",
  cancer2 = "Hematology/Oncology \\(R\\)/Cancer/",
  leukemia2 = "Hematology/Oncology \\(R\\)/Cancer/Hematologic Malignancy/(AML|ALL|CLL|CML|leukemia)",
  lymphoma2 = "Hematology/Oncology \\(R\\)/Cancer/Hematologic Malignancy/(non-Hodgkins lymphoma|Hodgkins disease)",
  mi1 = "Cardiovascular \\(R\\)/Myocardial Infarction/",
  chf1 = "Cardiovascular \\(R\\)/Congestive Heart Failure/",
  pvd1 = "Peripheral Vascular Disease/peripheral vascular disease$",
  tia1 = "Neurologic/TIA\\(s\\)/",
  dementia1 = "Neurologic/Dementia/dementia$",
  copd1 = "Pulmonary/COPD/",
  ctd1 = "Rheumatic/(SLE|Rheumatoid Arthritis|Scleroderma|Vasculitis|Dermato/Polymyositis)",
  pud1 = "Gastrointestinal \\(R\\)/Peptic Ulcer Disease/",
  liver1 = "Gastrointestinal \\(R\\)/Cirrhosis/(clinical diagnosis|biopsy proven)"
)

cci_weights <- c(
  mets6 = 6,
  aids6 = 6,
  liver3 = 3,
  stroke2 = 2,
  renal2 = 2,
  dm1 = 1,
  cancer2 = 2,
  leukemia2 = 2,
  lymphoma2 = 2,
  mi1 = 1,
  chf1 = 1,
  pvd1 = 1,
  tia1 = 1,
  dementia1 = 1,
  copd1 = 1,
  ctd1 = 1,
  pud1 = 1,
  liver1 = 1
)

ph <- pasthistory %>%
  semi_join(patient_cohort, by = "patientunitstayid") %>%
  mutate(
    pasthistorypath = as.character(pasthistorypath)
  )

if (nrow(ph) > 0) {
  
  ph_flags <- ph %>%
    transmute(
      patientunitstayid,
      mets6 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$mets6, ignore_case = TRUE))),
      aids6 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$aids6, ignore_case = TRUE))),
      liver3 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$liver3, ignore_case = TRUE))),
      stroke2 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$stroke2, ignore_case = TRUE))),
      renal2 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$renal2, ignore_case = TRUE))),
      dm1 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$dm1, ignore_case = TRUE))),
      cancer2 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$cancer2, ignore_case = TRUE))),
      leukemia2 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$leukemia2, ignore_case = TRUE))),
      lymphoma2 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$lymphoma2, ignore_case = TRUE))),
      mi1 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$mi1, ignore_case = TRUE))),
      chf1 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$chf1, ignore_case = TRUE))),
      pvd1 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$pvd1, ignore_case = TRUE))),
      tia1 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$tia1, ignore_case = TRUE))),
      dementia1 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$dementia1, ignore_case = TRUE))),
      copd1 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$copd1, ignore_case = TRUE))),
      ctd1 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$ctd1, ignore_case = TRUE))),
      pud1 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$pud1, ignore_case = TRUE))),
      liver1 = as.integer(str_detect(pasthistorypath, regex(cci_patterns$liver1, ignore_case = TRUE)))
    )
  
  cci_by_patient <- ph_flags %>%
    group_by(patientunitstayid) %>%
    summarise(
      across(
        all_of(names(cci_weights)),
        ~ max(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    mutate(
      charlson_score =
        mets6 * 6 +
        aids6 * 6 +
        liver3 * 3 +
        stroke2 * 2 +
        renal2 * 2 +
        dm1 * 1 +
        cancer2 * 2 +
        leukemia2 * 2 +
        lymphoma2 * 2 +
        mi1 * 1 +
        chf1 * 1 +
        pvd1 * 1 +
        tia1 * 1 +
        dementia1 * 1 +
        copd1 * 1 +
        ctd1 * 1 +
        pud1 * 1 +
        liver1 * 1
    ) %>%
    select(patientunitstayid, all_of(names(cci_weights)), charlson_score)
  
} else {
  
  cci_by_patient <- patient_cohort %>%
    select(patientunitstayid) %>%
    mutate(
      mets6 = 0L,
      aids6 = 0L,
      liver3 = 0L,
      stroke2 = 0L,
      renal2 = 0L,
      dm1 = 0L,
      cancer2 = 0L,
      leukemia2 = 0L,
      lymphoma2 = 0L,
      mi1 = 0L,
      chf1 = 0L,
      pvd1 = 0L,
      tia1 = 0L,
      dementia1 = 0L,
      copd1 = 0L,
      ctd1 = 0L,
      pud1 = 0L,
      liver1 = 0L,
      charlson_score = 0L
    )
}

age_adjust <- patient_cohort %>%
  select(patientunitstayid, age) %>%
  mutate(
    age_score = case_when(
      age >= 80 ~ 4L,
      age >= 70 ~ 3L,
      age >= 60 ~ 2L,
      age >= 50 ~ 1L,
      TRUE ~ 0L
    )
  )

charlson_df <- patient_cohort %>%
  select(patientunitstayid) %>%
  left_join(cci_by_patient, by = "patientunitstayid") %>%
  mutate(
    across(
      all_of(names(cci_weights)),
      ~ coalesce(.x, 0L)
    ),
    charlson_score = coalesce(charlson_score, 0L)
  ) %>%
  left_join(age_adjust, by = "patientunitstayid") %>%
  mutate(
    charlson_age_adjusted_score = charlson_score + age_score
  )

if (USE_AGE_ADJUSTED_CHARLSON) {
  charlson_df <- charlson_df %>%
    mutate(final_charlson_score = charlson_age_adjusted_score)
} else {
  charlson_df <- charlson_df %>%
    mutate(final_charlson_score = charlson_score)
}

charlson_component_counts <- charlson_df %>%
  summarise(
    across(
      all_of(names(cci_weights)),
      ~ sum(.x, na.rm = TRUE)
    ),
    charlson_score_mean = mean(charlson_score, na.rm = TRUE),
    charlson_score_sd = sd(charlson_score, na.rm = TRUE),
    charlson_score_median = median(charlson_score, na.rm = TRUE),
    charlson_score_q1 = quantile(charlson_score, 0.25, na.rm = TRUE),
    charlson_score_q3 = quantile(charlson_score, 0.75, na.rm = TRUE),
    charlson_score_max = max(charlson_score, na.rm = TRUE)
  )

############################################################
# 15b. Charlson/pastHistory audit
############################################################

if (RUN_CHARLSON_AUDIT) {
  ph_basic <- tibble(
    final_cohort_before_complete_case_n = nrow(patient_cohort),
    patients_with_any_pastHistory = n_distinct(ph$patientunitstayid),
    pastHistory_rows = nrow(ph),
    unique_pastHistory_paths = n_distinct(ph$pasthistorypath)
  )
  
  unique_paths <- ph %>%
    count(pasthistorypath, sort = TRUE)
  
  write_csv(unique_paths, OUT_PASTHISTORY_PATHS)
  
  pattern_hits <- purrr::imap_dfr(
    cci_patterns,
    function(pattern, name) {
      hit_idx <- str_detect(ph$pasthistorypath, regex(pattern, ignore_case = TRUE))
      
      tibble(
        charlson_component = name,
        pattern = pattern,
        matching_rows = sum(hit_idx, na.rm = TRUE),
        matching_patients = n_distinct(ph$patientunitstayid[hit_idx])
      )
    }
  )
  
  broad_terms <- c(
    "metast", "aids", "hiv", "cirrhos", "liver", "stroke", "cva",
    "renal", "kidney", "diabetes", "cancer", "leuk", "lymph",
    "myocardial", "infarct", "heart failure", "chf",
    "vascular", "dementia", "copd", "ulcer", "rheum"
  )
  
  broad_hits <- purrr::map_dfr(
    broad_terms,
    function(term) {
      hits <- ph %>%
        filter(str_detect(pasthistorypath, regex(term, ignore_case = TRUE)))
      
      tibble(
        term = term,
        matching_rows = nrow(hits),
        matching_patients = n_distinct(hits$patientunitstayid)
      )
    }
  )
  
  relevant_examples <- ph %>%
    filter(
      str_detect(
        pasthistorypath,
        regex(paste(broad_terms, collapse = "|"), ignore_case = TRUE)
      )
    ) %>%
    count(pasthistorypath, sort = TRUE) %>%
    head(100)
  
  sink(OUT_CHARLSON_AUDIT)
  
  cat("\n====================\n")
  cat("CHARLSON AUDIT V2\n")
  cat("====================\n\n")
  
  cat("====================\n")
  cat("BASIC PASTHISTORY COVERAGE\n")
  cat("====================\n")
  print(ph_basic)
  
  cat("\n====================\n")
  cat("CURRENT CHARLSON PATTERN HITS\n")
  cat("====================\n")
  print(pattern_hits, n = nrow(pattern_hits))
  
  cat("\n====================\n")
  cat("BROAD KEYWORD HITS\n")
  cat("====================\n")
  print(broad_hits, n = nrow(broad_hits))
  
  cat("\n====================\n")
  cat("TOP 100 RELEVANT PATH EXAMPLES\n")
  cat("====================\n")
  print(relevant_examples, n = 100)
  
  cat("\n====================\n")
  cat("TOP 100 ALL PASTHISTORY PATHS\n")
  cat("====================\n")
  print(head(unique_paths, 100), n = 100)
  
  sink()
}

############################################################
# Join Charlson scores to patient_cohort
# Remove any previous Charlson columns before joining
############################################################

patient_cohort <- patient_cohort %>%
  select(
    -any_of(c(
      names(cci_weights),
      "charlson_score",
      "charlson_age_adjusted_score",
      "final_charlson_score",
      "age_score"
    ))
  ) %>%
  left_join(
    charlson_df %>%
      select(
        patientunitstayid,
        charlson_score,
        charlson_age_adjusted_score,
        final_charlson_score
      ),
    by = "patientunitstayid"
  )

############################################################
# Safety check: Charlson should not be all zero
############################################################

if (all(patient_cohort$final_charlson_score == 0, na.rm = TRUE)) {
  stop("ERROR: final_charlson_score is all zero after joining Charlson scores.")
}

############################################################
# 16. Ethnicity harmonization
############################################################

patient_cohort <- patient_cohort %>%
  mutate(
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
    )
  )

############################################################
# 17. Final analytic dataset with explicit flowchart steps
############################################################

p7 <- patient_cohort %>%
  filter(!is.na(apacheiv))

add_flow("Available APACHE IV", p7)

p8 <- p7 %>%
  filter(!is.na(icu_mortality))

add_flow("Available ICU mortality", p8)

final_df <- p8 %>%
  mutate(
    icu_los_days = unitdischargeoffset / 1440
  ) %>%
  select(
    patientunitstayid,
    icu_mortality,
    er_los_hrs,
    vasopressors,
    age,
    sex,
    ethnicity,
    apacheadmissiondx,
    apacheiv,
    charlson_score,
    charlson_age_adjusted_score,
    final_charlson_score,
    icu_los_days
  ) %>%
  filter(
    !is.na(icu_mortality),
    !is.na(er_los_hrs),
    !is.na(age),
    !is.na(sex),
    !is.na(ethnicity),
    !is.na(vasopressors),
    !is.na(apacheiv),
    !is.na(final_charlson_score)
  )

add_flow("Final analytic cohort", final_df)

############################################################
# 18. Final checks
############################################################

flow <- flow %>%
  mutate(
    excluded_from_previous = lag(N) - N
  )

missing_summary <- final_df %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "n_missing"
  )

diagnosis_counts <- final_df %>%
  count(apacheadmissiondx, sort = TRUE)

ethnicity_counts <- final_df %>%
  count(ethnicity, icu_mortality, name = "n") %>%
  arrange(ethnicity, icu_mortality)

vaso_mortality <- final_df %>%
  count(vasopressors, icu_mortality, name = "n") %>%
  arrange(vasopressors, icu_mortality)

mortality_counts <- final_df %>%
  count(icu_mortality, name = "n") %>%
  mutate(percent = 100 * n / sum(n))

ed_los_summary <- final_df %>%
  summarise(
    n = n(),
    mean = mean(er_los_hrs),
    sd = sd(er_los_hrs),
    median = median(er_los_hrs),
    q1 = quantile(er_los_hrs, 0.25),
    q3 = quantile(er_los_hrs, 0.75),
    min = min(er_los_hrs),
    max = max(er_los_hrs)
  )

############################################################
# 19. Save outputs
############################################################

write_csv(final_df, OUT_CSV_COHORT)
write_csv(flow, OUT_FLOWCHART)
write_csv(diagnosis_counts, OUT_DX_COUNTS)
write_csv(missing_summary, OUT_MISSING)

############################################################
# 20. Console/audit report
############################################################

sink(OUT_AUDIT_TXT)

cat("\n====================\n")
cat("AUDIT SUMMARY V2\n")
cat("====================\n\n")

cat("Options:\n")
cat("USE_AGE_ADJUSTED_CHARLSON =", USE_AGE_ADJUSTED_CHARLSON, "\n")
cat("FILTER_UNITVISITNUMBER_1 =", FILTER_UNITVISITNUMBER_1, "\n")
cat("RUN_CHARLSON_AUDIT =", RUN_CHARLSON_AUDIT, "\n\n")

cat("====================\n")
cat("FLOWCHART\n")
cat("====================\n")
print(flow)

cat("\n====================\n")
cat("OFFSET AUDIT BEFORE ED LOS FILTER\n")
cat("====================\n")
print(offset_audit_before)

cat("\n====================\n")
cat("FINAL N\n")
cat("====================\n")
print(nrow(final_df))

cat("\n====================\n")
cat("ICU MORTALITY\n")
cat("====================\n")
print(mortality_counts)

cat("\n====================\n")
cat("VASOPRESSORS BY ICU MORTALITY\n")
cat("====================\n")
print(vaso_mortality)

cat("\n====================\n")
cat("ETHNICITY BY ICU MORTALITY\n")
cat("====================\n")
print(ethnicity_counts)

cat("\n====================\n")
cat("ED LOS SUMMARY\n")
cat("====================\n")
print(ed_los_summary)

cat("\n====================\n")
cat("MISSING VALUES IN FINAL DATASET\n")
cat("====================\n")
print(missing_summary)

cat("\n====================\n")
cat("CHARLSON COMPONENT COUNTS\n")
cat("====================\n")
print(charlson_component_counts)

cat("\n====================\n")
cat("CHARLSON SCORE SUMMARY IN CHARLSON_DF\n")
cat("====================\n")
print(summary(charlson_df$final_charlson_score))

cat("\n====================\n")
cat("CHARLSON SCORE SUMMARY IN PATIENT_COHORT AFTER JOIN\n")
cat("====================\n")
print(summary(patient_cohort$final_charlson_score))

cat("\n====================\n")
cat("CHARLSON SCORE SUMMARY IN FINAL DATASET\n")
cat("====================\n")
print(summary(final_df$final_charlson_score))

cat("\n====================\n")
cat("APACHE DUPLICATE AUDIT\n")
cat("====================\n")
print(apache_dup_audit)

cat("\n====================\n")
cat("APACHE VERSION AUDIT\n")
cat("====================\n")
print(apache_version_audit)

if ("unitvisitnumber" %in% names(patient_raw)) {
  cat("\n====================\n")
  cat("UNITVISITNUMBER AUDIT IN PATIENT_RAW\n")
  cat("====================\n")
  print(table(patient_raw$unitvisitnumber, useNA = "ifany"))
}

cat("\n====================\n")
cat("TOP 50 DIAGNOSES\n")
cat("====================\n")
print(head(diagnosis_counts, 50))

sink()

cat("\n====================\n")
cat("SCRIPT FINISHED\n")
cat("====================\n")
cat("Final dataset saved to:\n", OUT_CSV_COHORT, "\n")
cat("Flowchart saved to:\n", OUT_FLOWCHART, "\n")
cat("Diagnosis counts saved to:\n", OUT_DX_COUNTS, "\n")
cat("Missing summary saved to:\n", OUT_MISSING, "\n")
cat("Audit summary saved to:\n", OUT_AUDIT_TXT, "\n")
if (RUN_CHARLSON_AUDIT) {
  cat("Charlson audit saved to:\n", OUT_CHARLSON_AUDIT, "\n")
  cat("PastHistory unique paths saved to:\n", OUT_PASTHISTORY_PATHS, "\n")
}
cat("\n")

