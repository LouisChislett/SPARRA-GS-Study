# =============================================================================
# exclusions.R
#
# Author       : Louis Chislett
#
# Purpose :
#   Run missingness and exclusions-focused quality checks on the monthly
#   covariate matrix with target and write summary outputs for review.
#
#   IMPORTANT IMPLEMENTATION DETAIL:
#     - This script assumes covariate_matrix_monthly_with_target.fst already
#       exists and contains year_month, target, age, decile, sexM,
#       and subcohort_primary.
# =============================================================================

# Load libraries
library(fst)
library(dplyr)
library(purrr)
library(tibble)

# Read monthly modelling dataset
dir_clean <- "PATH/cleanData"
in_file   <- file.path(dir_clean, "covariate_matrix_monthly_with_target.fst")
df_all    <- read_fst(in_file)

# -----------------------------------------------------------------------------
# 1) Summary: counts and % TRUE/NA per month - sanity checking
# -----------------------------------------------------------------------------
summary_monthly <- df_all %>%
  group_by(year_month) %>%
  summarise(
    total       = n(),
    n_NA        = sum(is.na(target)),
    positives   = sum(target, na.rm = TRUE),
    negatives   = sum(!target, na.rm = TRUE),
    pct_true    = 100 * positives / (positives + negatives),
    pct_NA      = 100 * n_NA / total,
    .groups     = "drop"
  ) %>%
  arrange(year_month)

print(head(summary_monthly, 12))   # first year
print(tail(summary_monthly, 12))   # last year


# =============================================================================
# 2) Missingness / age<30 combination checks for SPECIFIC months:
#     Jan-2012 and Oct-2022
# =============================================================================

required_cols <- c("year_month", "target", "age", "decile", "sexM")
missing_cols <- setdiff(required_cols, names(df_all))
if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Cannot run combination checks; missing columns in df_all: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# flags: 4 missingness flags + 1 age-under-30 flag
df_flags <- df_all %>%
  transmute(
    year_month,
    NA_target   = is.na(target),
    NA_age      = is.na(age),
    NA_decile   = is.na(decile),
    NA_sexM     = is.na(sexM),
    age_under30 = !is.na(age) & age < 30
  )

flag_names <- c("NA_target", "NA_decile", "NA_age", "NA_sexM", "age_under30")

# helper: count rows where ALL selected flags are TRUE
count_all_true <- function(df, selected_flag_names) {
  if (length(selected_flag_names) == 0) {
    return(nrow(df))
  }
  mat <- as.data.frame(df[, selected_flag_names, drop = FALSE])
  sum(rowSums(mat) == length(selected_flag_names))
}

# helper: make output labels look nice
make_combination_label <- function(sel) {
  if (length(sel) == 0) {
    return("none_missing")
  }
  
  sel_clean <- gsub("^NA_target$", "died", sel)
  sel_clean <- gsub("^NA_age$", "age", sel_clean)
  sel_clean <- gsub("^NA_decile$", "decile", sel_clean)
  sel_clean <- gsub("^NA_sexM$", "sex", sel_clean)
  sel_clean <- gsub("^age_under30$", "age_lt_30", sel_clean)
  
  paste(sel_clean, collapse = "_and_")
}

# Build summary for a single month label (YYYY-MM)
summarise_one_month <- function(month_label) {
  dfm <- df_flags %>% filter(year_month == month_label)
  
  if (nrow(dfm) == 0) {
    return(
      tibble(
        year_month  = month_label,
        combination = character(),
        n           = integer(),
        pct         = numeric()
      )
    )
  }
  
  out <- map_dfr(0:length(flag_names), function(k) {
    combos <- if (k == 0) list(character(0)) else combn(flag_names, k, simplify = FALSE)
    
    map_dfr(combos, function(sel) {
      n <- count_all_true(dfm, sel)
      
      tibble(
        year_month  = month_label,
        combination = make_combination_label(sel),
        n           = n,
        pct         = 100 * n / nrow(dfm)
      )
    })
  }) %>%
    arrange(combination)
  
  out
}

jan2012 <- summarise_one_month("2012-01")
oct2022 <- summarise_one_month("2022-10")

cat("\n--- Missingness / age<30 combination checks: Jan-2012 ---\n")
print(jan2012)

cat("\n--- Missingness / age<30 combination checks: Oct-2022 ---\n")
print(oct2022)

# write out the two tables to ./exclusions
dir_exclusions <- file.path(getwd(), "exclusions")
if (!dir.exists(dir_exclusions)) dir.create(dir_exclusions, recursive = TRUE)

write.csv(jan2012, file.path(dir_exclusions, "combinations_2012-01.csv"), row.names = FALSE)
write.csv(oct2022, file.path(dir_exclusions, "combinations_2022-10.csv"), row.names = FALSE)

# =============================================================================
# 3) Counts of eligible people with each subcohort indicator
#     after exclusions, for 2012-01 and 2022-10
# =============================================================================

required_cols_3 <- c(
  "year_month", "target", "age", "decile", "sexM",
  "subcohort_FE", "subcohort_LTC", "subcohort_YED"
)

missing_cols_3 <- setdiff(required_cols_3, names(df_all))
if (length(missing_cols_3) > 0) {
  stop(
    paste0(
      "Cannot run eligible subcohort summaries; missing columns in df_all: ",
      paste(missing_cols_3, collapse = ", ")
    )
  )
}

months_to_check <- c("2012-01", "2022-10")

# eligible = survives the exclusion criteria
df_eligible <- df_all %>%
  filter(
    !is.na(target),   # exclude died / missing target
    !is.na(age),
    age >= 30,
    !is.na(decile),
    !is.na(sexM)
  )

summarise_subcohort_indicators <- function(data, month_label) {
  dfm <- data %>%
    filter(year_month == month_label)
  
  if (nrow(dfm) == 0) {
    return(tibble(
      year_month = character(),
      indicator  = character(),
      n          = integer(),
      prop       = numeric(),
      pct        = numeric()
    ))
  }
  
  tibble(
    year_month = month_label,
    indicator  = c("subcohort_FE", "subcohort_LTC", "subcohort_YED"),
    n = c(
      sum(dplyr::coalesce(dfm$subcohort_FE, 0) == 1),
      sum(dplyr::coalesce(dfm$subcohort_LTC, 0) == 1),
      sum(dplyr::coalesce(dfm$subcohort_YED, 0) == 1)
    )
  ) %>%
    mutate(
      prop = n / nrow(dfm),
      pct  = 100 * prop
    )
}

subcohort_indicator_summary <- map_dfr(
  months_to_check,
  ~ summarise_subcohort_indicators(df_eligible, .x)
)

cat("\n--- Eligible people with each subcohort indicator: 2012-01 and 2022-10 ---\n")
print(subcohort_indicator_summary)

write.csv(
  subcohort_indicator_summary,
  file.path(dir_exclusions, "eligible_subcohort_indicator_counts_2012-01_2022-10.csv"),
  row.names = FALSE
)

# =============================================================================
# 4) Overlaps between eligible subcohort indicators
#     after exclusions, for 2012-01 and 2022-10
# =============================================================================

summarise_subcohort_overlaps <- function(data, month_label) {
  dfm <- data %>%
    filter(year_month == month_label) %>%
    mutate(
      FE  = dplyr::coalesce(subcohort_FE, 0)  == 1,
      LTC = dplyr::coalesce(subcohort_LTC, 0) == 1,
      YED = dplyr::coalesce(subcohort_YED, 0) == 1
    )
  
  if (nrow(dfm) == 0) {
    return(tibble(
      year_month   = character(),
      combination  = character(),
      n            = integer(),
      prop         = numeric(),
      pct          = numeric()
    ))
  }
  
  total_n <- nrow(dfm)
  
  tibble(
    year_month  = month_label,
    combination = c(
      "total_unique_individuals",
      "none",
      "FE_only",
      "LTC_only",
      "YED_only",
      "FE_and_LTC_only",
      "FE_and_YED_only",
      "LTC_and_YED_only",
      "FE_and_LTC_and_YED"
    ),
    n = c(
      total_n,
      sum(!dfm$FE & !dfm$LTC & !dfm$YED),
      sum( dfm$FE & !dfm$LTC & !dfm$YED),
      sum(!dfm$FE &  dfm$LTC & !dfm$YED),
      sum(!dfm$FE & !dfm$LTC &  dfm$YED),
      sum( dfm$FE &  dfm$LTC & !dfm$YED),
      sum( dfm$FE & !dfm$LTC &  dfm$YED),
      sum(!dfm$FE &  dfm$LTC &  dfm$YED),
      sum( dfm$FE &  dfm$LTC &  dfm$YED)
    )
  ) %>%
    mutate(
      prop = n / total_n,
      pct  = 100 * prop
    )
}

subcohort_overlap_summary <- map_dfr(
  months_to_check,
  ~ summarise_subcohort_overlaps(df_eligible, .x)
)

cat("\n--- Overlaps between eligible subcohort indicators: 2012-01 and 2022-10 ---\n")
print(subcohort_overlap_summary)

write.csv(
  subcohort_overlap_summary,
  file.path(dir_exclusions, "eligible_subcohort_overlap_counts_2012-01_2022-10.csv"),
  row.names = FALSE
)