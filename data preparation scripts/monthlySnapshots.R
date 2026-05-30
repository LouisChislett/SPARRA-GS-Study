# =============================================================================
# monthlySnapshots.R
#
# Author       : Louis Chislett
#
# Purpose :
#   Build a monthly SPARRA v3-like covariate matrix with a forward-looking
#   target and write it to covariate_matrix_monthly_with_target.fst.
#
#   This script:
#     (1) Reads cleaned patient-level and event-level source data
#         - Loads .fst inputs from the cleanData directory
#         - Bundles source tables for downstream feature and target generation
#
#     (2) Builds monthly patient snapshots across the study window
#         - Generates monthly features with transformer_v3like()
#         - Expands to a complete patient × month skeleton for each cutoff
#
#     (3) Applies feature filling rules and derives the outcome target
#         - Fills missing event-derived features with 0/FALSE after the join
#         - Preserves missing demographic values for age, decile, and sexM
#         - Computes the 12-month target using add_target()
#
#     (4) Combines and writes the final modelling dataset
#         - Row-binds all monthly snapshots into a single dataset
#         - Writes covariate_matrix_monthly_with_target.fst to disk
#
#   IMPORTANT IMPLEMENTATION DETAIL:
#     - Missing target values are retained; in this pipeline they indicate
#       patients who died before the study period.
#     - Missingness checks, exclusion summaries, and downstream QC are handled
#       separately in exclusions.R.
# =============================================================================

# Load libraries
library(fst)
library(dplyr)
library(lubridate)
library(purrr)
library(tibble)
library(tidyr)

# Source helper functions
source("transformerV3Like.R")   # feature transformer (your updated version)
source("addTarget.R")           # add_target()
source("groupings.R")           # grouping utils
source("icd10ChaptersLookup.R")


# -----------------------------------------------------------------------------
# Read cleaned data
# -----------------------------------------------------------------------------
dir_clean <- "PATH/cleanData"
patients  <- read_fst(file.path(dir_clean, "patients.fst"))
episodes  <- read_fst(file.path(dir_clean, "episodes.fst"))
AE2       <- read_fst(file.path(dir_clean, "AE2.fst"))
deaths    <- read_fst(file.path(dir_clean, "deaths.fst"))
PIS       <- read_fst(file.path(dir_clean, "PIS.fst"))
SMR00     <- read_fst(file.path(dir_clean, "SMR00.fst"))
SMR01     <- read_fst(file.path(dir_clean, "SMR01.fst"))
SMR04     <- read_fst(file.path(dir_clean, "SMR04.fst"))
SPARRALTC <- read_fst(file.path(dir_clean, "SPARRALTC.fst"))

# Bundle data tables
list_of_data_tables <- list(
  AE2       = AE2,
  deaths    = deaths,
  PIS       = PIS,
  SMR00     = SMR00,
  SMR01     = SMR01,
  SMR04     = SMR04,
  SPARRALTC = SPARRALTC,
  episodes  = episodes
)

# -----------------------------------------------------------------------------
# Variables to keep from transformer output
# IMPORTANT: include age_missing / sex_missing / decile_missing if you want them
# to persist in the output dataset (df_all).
# -----------------------------------------------------------------------------
v3_vars <- c(
  "id","time","target",
  "subcohort_FE","subcohort_LTC","subcohort_YED","subcohort_primary",
  "num_emergency_admissions","emergency_bed_days",
  "num_elective_admissions","num_daycase_admissions","elective_bed_days",
  "num_alcohol_drug_admissions","num_alcohol_admissions",
  "numLTCs_resulting_in_admin",
  "num_bnf_sections","pis_respiratory","pis_cns","pis_infections",
  "pis_endocrine","pis_incontinence","pis_sub_depend","pis_dementia",
  "pis_corticosteroids","pis_fluids","pis_nutrition","pis_vitamins",
  "pis_bandages","pis_catheters","pis_stoma","pis_gut_motility",
  "pis_antisecretory","pis_intestinal","pis_anticoagulant",
  "pis_antifibrinolytic","pis_mucolytics","pis_diabetes",
  "num_psych_admissions",
  "num_ae2_attendances",
  "num_outpatient_appointment_general","num_outpatient_appointment_psych",
  
  # demographics + missingness flags
  "age","age_missing",
  "decile","decile_missing",
  "sexM","sex_missing",
  
  "parkinsons_indicated","MS_indicated","epilepsy_indicated",
  "dementia_indicated",
  "first_event_reason","first_event_time","emerg_reason","emerg_first_time"
)

# -----------------------------------------------------------------------------
# Define monthly cutoffs: 2012-01-01 to 2022-10-01 (inclusive)
# -----------------------------------------------------------------------------
start_date   <- as_date("2012-01-01")
end_date     <- as_date("2022-10-01")
time_cutoffs <- as_datetime(seq(from = start_date, to = end_date, by = "month"), tz = "UTC")

# -----------------------------------------------------------------------------
# Iterate over each cutoff, compute features + target
# -----------------------------------------------------------------------------
datasets <- map(time_cutoffs, function(cutoff) {
  cat("Processing snapshot:", format(cutoff, "%Y-%m-%d"), "\n")
  
  # Generate raw features for this cutoff
  out_feats <- transformer_v3like(
    patients            = patients,
    episodes            = episodes,
    list_of_data_tables = list_of_data_tables,
    time_cutoffs        = as.numeric(cutoff),
    keep_vars           = v3_vars
  )
  
  feat_mat <- if (is.list(out_feats) && "matrix" %in% names(out_feats)) {
    out_feats$matrix
  } else {
    out_feats
  }
  
  # Convert to tibble, preserving `id`
  if (is.matrix(feat_mat)) {
    feat_df <- as_tibble(feat_mat, rownames = "id") %>%
      mutate(id = as.integer(id))
  } else {
    feat_df <- as_tibble(feat_mat)
  }
  
  # De-dupe per id (in case transformer returns repeats)
  feat_df <- feat_df %>% distinct(id, .keep_all = TRUE)
  
  # Build full patient × time skeleton
  base_df <- patients %>%
    select(id) %>%
    distinct() %>%
    mutate(
      time       = cutoff,
      year       = year(cutoff),
      month      = month(cutoff),
      year_month = format(cutoff, "%Y-%m")
    )
  
  # Left-join features onto skeleton
  full_feat_df <- base_df %>%
    left_join(feat_df, by = "id")
  
  # IMPORTANT FILL RULE:
  #   - Fill event-derived features with 0/FALSE
  #   - DO NOT overwrite missing demographics age/decile/sexM (keep NA)
  protected_demo <- c("age", "decile", "sexM")
  
  # Fill logical columns (safe to fill FALSE)
  full_feat_df <- full_feat_df %>%
    mutate(across(where(is.logical), ~ replace_na(.x, FALSE)))
  
  # Fill integer columns EXCEPT protected demographics
  int_cols <- names(full_feat_df)[vapply(full_feat_df, is.integer, logical(1))]
  int_fill <- setdiff(int_cols, protected_demo)
  if (length(int_fill) > 0) {
    full_feat_df <- full_feat_df %>%
      mutate(across(all_of(int_fill), ~ replace_na(.x, 0L)))
  }
  
  # Fill numeric columns EXCEPT protected demographics (in case they are numeric)
  num_cols <- names(full_feat_df)[vapply(full_feat_df, is.numeric, logical(1))]
  num_fill <- setdiff(num_cols, protected_demo)
  if (length(num_fill) > 0) {
    full_feat_df <- full_feat_df %>%
      mutate(across(all_of(num_fill), ~ replace_na(.x, 0)))
  }
  
  # Handle the factor label explicitly (if present)
  if ("subcohort_primary" %in% names(full_feat_df)) {
    sp <- as.character(full_feat_df$subcohort_primary)
    sp[is.na(sp)] <- "None"
    full_feat_df$subcohort_primary <- factor(sp, levels = c("None","LTC","YED","FE"))
  }
  
  # Compute 12-month target from this monthly snapshot
  tgt_df <- add_target(
    data                 = full_feat_df %>% select(id, time),
    list_of_data_tables  = list_of_data_tables,
    period               = 365,
    ae2_window           = 2,
    include_reason       = FALSE,
    include_emerg_reason = TRUE
  )
  
  # Join target (target may be NA for died pre-study; keep it as NA)
  full_feat_df <- full_feat_df %>%
    left_join(
      tgt_df %>%
        select(id, time, target,
               first_event_time, first_event_reason,
               emerg_first_time, emerg_reason),
      by = c("id", "time")
    )
  
  full_feat_df
})

# Combine all months
df_all <- bind_rows(datasets)

# -----------------------------------------------------------------------------
# Write to disk
# -----------------------------------------------------------------------------
out_file <- file.path(dir_clean, "covariate_matrix_monthly_with_target.fst")
write_fst(df_all, out_file)
cat("Wrote combined monthly dataset to", out_file, "\n")
