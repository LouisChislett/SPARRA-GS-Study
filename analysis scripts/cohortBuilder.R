# =============================================================================
# cohortBuilder.R
#
# Author       : Louis Chislett
#
# Purpose :
#   Build SPARRA CPM v3 monthly snapshot datasets, including exclusions and adding
#   post stratification weights. Script includes:
#
#     (1) DIRECT STANDARDISATION BY AGE to a reference month distribution
#         - overall
#         - FE-specific
#         - LTC-specific
#         - YED-specific
#
#     (2) DIRECT STANDARDISATION TO EXTERNAL POPULATION
#         - using Scotland mid-year population counts
#         - weighting by:
#             * age band (5-year bands from 30; with 85-89 and 90+ collapsed
#               to 85+ for population weights only)
#             * sex (cohort variable: sexM, coded 1=Male, 0=Female)
#             * SIMD quintile (derived from decile for population weights only)
#         - computed on the FULL dataset only
#         - produces:
#             w_pop_raw, w_pop_stab, w_pop_trim
#
#   IMPORTANT IMPLEMENTATION DETAIL:
#     - Cohort-specific age weights are computed using ONLY rows in that cohort,
#       but are ATTACHED back to ALL rows (as NA for non-members).
#
# Key functions for downstream scripts:
#   - add_age_standardisation_weights(...)
#   - add_population_standardisation_weights(...)
#   - add_all_age_weights_overall_and_subcohort(...)
#   - add_all_weights(...)
#   - make_sparra_v3_datasets(df_in) - creates monthly snapshot datasets
#   - make_sparra_v3_datasets_raw(df_in) - creates monthly snapshot datasets without variable binning/categorisation
# =============================================================================


# =============================================================================
# 0) Constants
# =============================================================================

# ---- minimum age included anywhere ----
MIN_AGE_INCLUDED <- 30

# ---- adult age bands used for all standardisation ----
# Using right = FALSE => intervals are [a, b), so:
#   30–34 == [30,35)
#   35–39 == [35,40)
#   ...
#   85–89 == [85,90)
#   90+   == [90,Inf)
AGE_BREAKS_STD <- c(30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, Inf)
AGE_LABELS_STD <- c(
  "30-34",
  "35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74",
  "75-79","80-84","85-89","90+"
)

# ---- external population file ----
POP_STANDARDISATION_CSV <- "scotland_midyear_pop_age_sex_simd_decile_2012_2022_no_persons.csv"


# =============================================================================
# 1) Helper utilities
# =============================================================================
to_ym <- function(x) {
  suppressPackageStartupMessages({ library(lubridate) })
  
  if (inherits(x, "Date"))   return(format(floor_date(x, "month"), "%Y-%m"))
  if (inherits(x, "POSIXt")) return(format(floor_date(as.Date(x), "month"), "%Y-%m"))
  
  if (is.numeric(x)) {
    xdt <- as.POSIXct(x, origin = "1970-01-01", tz = "UTC")
    return(format(floor_date(as.Date(xdt), "month"), "%Y-%m"))
  }
  
  if (is.character(x)) {
    ok_ym  <- grepl("^\\d{4}-\\d{2}$", x)
    ok_ymd <- grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
    if (all(ok_ym | ok_ymd | is.na(x))) {
      return(ifelse(ok_ym, x, format(floor_date(as.Date(x), "month"), "%Y-%m")))
    }
    suppressWarnings(return(format(floor_date(as.Date(paste0(x, "-01")), "month"), "%Y-%m")))
  }
  
  stop("time_col must be Date/POSIXt/numeric or 'YYYY-MM'/'YYYY-MM-DD'.")
}

# sexM is explicitly coded:
#   1 = Male   -> "Males"
#   0 = Female -> "Females"
#   it is encoded slightly differently in the population reference, so we define two helper functions here
standardise_sex_from_cohort <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  dplyr::case_when(
    !is.na(x_num) & x_num == 1 ~ "Males",
    !is.na(x_num) & x_num == 0 ~ "Females",
    TRUE ~ NA_character_
  )
}

standardise_sex_from_population <- function(x) {
  y <- trimws(as.character(x))
  y_low <- tolower(y)
  
  dplyr::case_when(
    y_low %in% c("male", "males", "m") ~ "Males",
    y_low %in% c("female", "females", "f") ~ "Females",
    TRUE ~ NA_character_
  )
}

parse_population_age_to_numeric <- function(x) {
  y <- trimws(as.character(x))
  y_low <- tolower(y)

  out <- rep(NA_real_, length(y))

  is_total <- y_low %in% c("total", "all ages", "all_age", "all")
  is_plus  <- grepl("^\\d+\\+$", y)
  is_num   <- grepl("^\\d+$", y)
  is_range <- grepl("^\\d+\\s*[-–]\\s*\\d+$", y)

  out[is_plus]  <- suppressWarnings(as.numeric(sub("\\+$", "", y[is_plus])))
  out[is_num]   <- suppressWarnings(as.numeric(y[is_num]))
  out[is_range] <- suppressWarnings(as.numeric(sub("\\s*[-–].*$", "", y[is_range])))

  # Keep totals and anything else unparseable as NA without emitting warnings
  out[is_total] <- NA_real_

  out
}

# For POPULATION weighting only:
#   - collapse 85-89 and 90+ into a single 85+ band
#   - convert SIMD deciles to quintiles: (1,2)->1, (3,4)->2, ..., (9,10)->5
collapse_age_band_for_population_weights <- function(x) {
  x_chr <- as.character(x)
  x_chr[!is.na(x_chr) & x_chr %in% c("85-89", "90+")] <- "85+"
  factor(
    x_chr,
    levels = c(
      "30-34",
      "35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74",
      "75-79","80-84","85+"
    ),
    ordered = TRUE
  )
}

decile_to_quintile <- function(x) {
  x_num <- suppressWarnings(as.integer(x))
  dplyr::case_when(
    !is.na(x_num) & x_num >= 1 & x_num <= 10 ~ ((x_num - 1L) %/% 2L) + 1L,
    TRUE ~ NA_integer_
  )
}


# =============================================================================
# Generic direct standardisation by age
#   - post-stratification weights to reference month age distribution
#   - can be applied overall or within a cohort subset
#   - returns only trimmed + normalised weights: <prefix>_trim
#   - includes support check for missing reference age bands in later months
# =============================================================================
add_age_standardisation_weights <- function(
    df,
    time_col   = "time",
    age_col    = "age",
    ref_month  = "2012-01",
    subset_filter = NULL,      # optional dplyr expression evaluated in df
    prefix     = "w_age",      # creates <prefix>_trim
    age_breaks = AGE_BREAKS_STD,
    age_labels = AGE_LABELS_STD,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = 3.5
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
  })
  
  
  # ---- define output weight column and reference month ----
  # Example: prefix = "w_age_FE" creates "w_age_FE_trim".
  weight_col <- paste0(prefix, "_trim")
  ref_ym <- to_ym(ref_month)
  
  # ---- prepare analysis rows ----
  # Convert time to YYYY-MM, coerce age to numeric, create standard age bands,
  # and apply the global exclusion rule for missing age / age below threshold.
  df2 <- df %>%
    mutate(
      `..ym..` = to_ym(.data[[time_col]]),
      `..age_num..` = suppressWarnings(as.numeric(.data[[age_col]])),
      age_band_std = cut(
        `..age_num..`,
        breaks = age_breaks,
        labels = age_labels,
        right = FALSE,
        ordered_result = TRUE
      )
    ) %>%
    filter(
      !is.na(`..age_num..`),
      `..age_num..` >= min_age_included,
      !is.na(age_band_std)
    ) %>%
    mutate(
      age_band_std = factor(age_band_std, levels = age_labels, ordered = TRUE)
    )
  
  # ---- choose weighting source ----
  # If subset_filter is NULL, weights are calculated using the full dataset.
  # If subset_filter is supplied, weights are calculated only within that cohort
  # subset, e.g. FE-only, LTC-only, or YED-only.
  df_wsrc <- df2
  if (!is.null(subset_filter)) {
    df_wsrc <- df2 %>% filter(!!subset_filter)
  }
  
  # ---- handle empty weighting source ----
  # If no rows remain after exclusions/subsetting, return the filtered data with
  # the requested weight column set to NA.
  if (nrow(df_wsrc) == 0) {
    out <- df2 %>%
      mutate(!!weight_col := NA_real_) %>%
      select(-`..ym..`, -`..age_num..`)
    
    return(list(
      data = out,
      weights_table = tibble(),
      qa = tibble()
    ))
  }
  
  # ---- calculate reference-month age distribution ----
  # p_ref is the target age distribution:
  #   p_ref = number in age band during reference month /
  #           total number in reference month.
  #
  # For cohort-specific calls, this is the reference distribution within that
  # subcohort only.
  ref_tab <- df_wsrc %>%
    filter(`..ym..` == ref_ym) %>%
    count(age_band_std, name = "n_ref") %>%
    mutate(p_ref = n_ref / sum(n_ref)) %>%
    select(age_band_std, n_ref, p_ref)
  
  # ---- calculate observed monthly age distributions ----
  # p_m is the observed age distribution in each month:
  #   p_m = number in age band during month m /
  #         total number in month m.
  month_tab <- df_wsrc %>%
    count(`..ym..`, age_band_std, name = "n_m") %>%
    group_by(`..ym..`) %>%
    mutate(p_m = n_m / sum(n_m)) %>%
    ungroup()
  
  # ---------------------------------------------------------------------------
  # Support check:
  # Each month must contain every age band that has positive share in the
  # reference month. Otherwise p_ref / p_m would imply infinite weights.
  #
  # This is done by creating the full month x reference-age-band grid, joining
  # the observed monthly counts, and checking whether any expected age band is
  # absent in any month.
  # ---------------------------------------------------------------------------
  support_check <- df_wsrc %>%
    distinct(`..ym..`) %>%
    tidyr::crossing(
      ref_tab %>%
        filter(!is.na(p_ref), p_ref > 0) %>%
        select(age_band_std, p_ref)
    ) %>%
    left_join(
      month_tab %>% select(`..ym..`, age_band_std, n_m, p_m),
      by = c("..ym..", "age_band_std")
    )
  
  
  # ---- calculate raw post-stratification weights ----
  # Raw weight for month m and age band a:
  #   w_raw[m, a] = p_ref[a] / p_m[m, a]
  #
  # Age bands under-represented relative to the reference month get weight > 1.
  # Age bands over-represented relative to the reference month get weight < 1.
  wtab <- month_tab %>%
    left_join(
      ref_tab %>% select(age_band_std, p_ref),
      by = "age_band_std"
    ) %>%
    mutate(w_raw = p_ref / p_m)
  
  # ---- attach raw weights to individual rows ----
  # Each row receives the raw weight for its month and age band.
  out <- df2 %>%
    left_join(
      wtab %>% select(`..ym..`, age_band_std, w_raw),
      by = c("..ym..", "age_band_std")
    )
  
  # ---- blank cohort-specific weights for non-members ----
  # For cohort-specific calls, the weights are calculated within the cohort,
  # but the data returned still contains all eligible rows. Non-members get NA.
  if (!is.null(subset_filter)) {
    in_subset <- df2 %>%
      transmute(`..in_subset..` = !!subset_filter) %>%
      pull(`..in_subset..`)
    
    out$w_raw[!in_subset] <- NA_real_
  }
  
  # ---- trim and normalise weights ----
  # Step 1: compute a month-specific cap:
  #   cap[m] = weight_cap x median raw weight in month m.
  #
  # Step 2: trim raw weights at that cap:
  #   w_trim0 = min(w_raw, cap[m]).
  #
  # Step 3: normalise the trimmed weights so their mean is 1 within each month:
  #   w_trim = w_trim0 / mean(w_trim0 in month m).
  out <- out %>%
    group_by(`..ym..`) %>%
    mutate(
      `..cap..` = weight_cap * median(w_raw, na.rm = TRUE),
      `..w_trim..` = ifelse(
        is.na(w_raw),
        NA_real_,
        pmin(w_raw, `..cap..`)
      ),
      `..w_trim..` = ifelse(
        is.na(`..w_trim..`),
        NA_real_,
        `..w_trim..` / mean(`..w_trim..`, na.rm = TRUE)
      )
    ) %>%
    ungroup() %>%
    mutate(
      !!weight_col := `..w_trim..`
    )
  
  # ---- quality-assurance summary ----
  # For each month, report:
  #   n     = number of rows with non-missing final weight
  #   max_w = maximum final weight
  #   ess   = effective sample size based on final weights
  qa <- out %>%
    filter(!is.na(.data[[weight_col]])) %>%
    group_by(`..ym..`) %>%
    summarise(
      n = n(),
      max_w = max(.data[[weight_col]], na.rm = TRUE),
      ess = (sum(.data[[weight_col]], na.rm = TRUE)^2) /
        sum(.data[[weight_col]]^2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(`..ym..`)
  
  # ---- weight lookup table ----
  # This is useful for checking the underlying month x age-band proportions
  # and the raw weights before trimming/normalisation.
  weights_table <- wtab %>%
    select(`..ym..`, age_band_std, n_m, p_m, p_ref, w_raw) %>%
    arrange(`..ym..`, age_band_std)
  
  # ---- return cleaned data and diagnostics ----
  # Temporary calculation columns are removed from the returned data.
  list(
    data = out %>%
      select(
        -`..ym..`,
        -`..age_num..`,
        -w_raw,
        -`..cap..`,
        -`..w_trim..`
      ),
    weights_table = weights_table,
    qa = qa
  )
}


# =============================================================================
# External population standardisation by age band x sex x SIMD quintile
#   - uses year-specific Scotland mid-year population totals
#   - for any month YYYY-MM, weights are matched to population year YYYY
#   - computes full-dataset population weights only
#   - returns only trimmed + normalised weights: <prefix>_trim
# =============================================================================
add_population_standardisation_weights <- function(
    df,
    pop_csv    = POP_STANDARDISATION_CSV,
    time_col   = "time",
    age_col    = "age",
    sex_col    = "sexM",
    decile_col = "decile",
    prefix     = "w_pop",
    age_breaks = AGE_BREAKS_STD,
    age_labels = AGE_LABELS_STD,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = 3.5
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
  })
  
  # ---- define output weight column ----
  # Example: prefix = "w_pop" creates "w_pop_trim".
  weight_col <- paste0(prefix, "_trim")
  
  
  # ---- read external population data ----
  pop_raw <- suppressMessages(
    readr::read_csv(pop_csv, show_col_types = FALSE)
  )

  
  # ---- prepare cohort rows ----
  # Each cohort row is assigned:
  #   - month: YYYY-MM
  #   - year: YYYY, used to match to same-year mid-year population
  #   - numeric age
  #   - standardised sex label
  #   - SIMD quintile derived from decile
  #   - standard age band
  #
  # Rows with missing/ineligible age, sex, or SIMD are excluded from population
  # weighting because they cannot be assigned to a population stratum.
  df2 <- df %>%
    mutate(
      `..ym..` = to_ym(.data[[time_col]]),
      `..year..` = suppressWarnings(as.integer(substr(`..ym..`, 1, 4))),
      `..age_num..` = suppressWarnings(as.numeric(.data[[age_col]])),
      `..sex_std..` = standardise_sex_from_cohort(.data[[sex_col]]),
      `..decile_std..` = decile_to_quintile(.data[[decile_col]]),
      age_band_std = cut(
        `..age_num..`,
        breaks = age_breaks,
        labels = age_labels,
        right = FALSE,
        ordered_result = TRUE
      )
    ) %>%
    filter(
      !is.na(`..ym..`),
      !is.na(`..year..`),
      !is.na(`..age_num..`),
      `..age_num..` >= min_age_included,
      !is.na(age_band_std),
      !is.na(`..sex_std..`),
      !is.na(`..decile_std..`)
    ) %>%
    mutate(
      # For population weighting only, 85-89 and 90+ are collapsed to 85+.
      age_band_std = collapse_age_band_for_population_weights(age_band_std)
    )
  
  # ---- handle empty cohort data after exclusions ----
  if (nrow(df2) == 0) {
    out <- df2 %>%
      mutate(!!weight_col := NA_real_) %>%
      select(
        -`..ym..`,
        -`..year..`,
        -`..age_num..`,
        -`..sex_std..`,
        -`..decile_std..`
      )
    
    return(list(
      data = out,
      weights_table = tibble(),
      qa = tibble()
    ))
  }
  
  # ---- prepare external population table ----
  # Convert population records to the same strata used in the cohort:
  #   year x age band x sex x SIMD quintile.
  #
  # p_ref is the external population share in each stratum within each year:
  #   p_ref = population count in stratum / total population count in year.
  pop_tab <- pop_raw %>%
    mutate(
      year = suppressWarnings(as.integer(year)),
      sex_std = standardise_sex_from_population(sex),
      simd_decile = suppressWarnings(as.integer(simd_decile)),
      simd_quintile = decile_to_quintile(simd_decile),
      age_num = parse_population_age_to_numeric(age),
      population = suppressWarnings(as.numeric(population))
    ) %>%
    filter(
      !is.na(year),
      !is.na(sex_std),
      !is.na(simd_quintile),
      !is.na(age_num),
      age_num >= min_age_included,
      !is.na(population)
    ) %>%
    mutate(
      age_band_std = cut(
        age_num,
        breaks = age_breaks,
        labels = age_labels,
        right = FALSE,
        ordered_result = TRUE
      ),
      age_band_std = collapse_age_band_for_population_weights(age_band_std)
    ) %>%
    filter(!is.na(age_band_std)) %>%
    group_by(year, age_band_std, sex_std, simd_quintile) %>%
    summarise(
      pop_n = sum(population, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(year) %>%
    mutate(
      p_ref = pop_n / sum(pop_n)
    ) %>%
    ungroup()
  
  # ---- calculate observed monthly cohort distribution over joint strata ----
  # For each month, calculate the cohort share in each:
  #   age band x sex x SIMD quintile
  #
  # p_m is the observed cohort share in that month:
  #   p_m = cohort count in stratum during month m /
  #         total cohort count in month m.
  month_tab <- df2 %>%
    count(
      `..ym..`,
      `..year..`,
      age_band_std,
      `..sex_std..`,
      `..decile_std..`,
      name = "n_m"
    ) %>%
    group_by(`..ym..`) %>%
    mutate(
      p_m = n_m / sum(n_m)
    ) %>%
    ungroup()
  
  # ---------------------------------------------------------------------------
  # Support check:
  # For each month, every population stratum with positive population share in
  # that year must also exist in the cohort month. Otherwise p_ref / p_m would
  # require division by zero and imply infinite weights.
  #
  # This creates all expected month x population-stratum combinations for the
  # relevant year, then joins the observed cohort month-level counts.
  # ---------------------------------------------------------------------------
  month_years <- df2 %>%
    distinct(`..ym..`, `..year..`)
  
  pop_tab_join <- pop_tab %>%
    rename(
      `..year..` = year,
      `..sex_std..` = sex_std,
      `..decile_std..` = simd_quintile
    )
  
  full_grid <- merge(
    month_years,
    pop_tab_join,
    by = "..year.."
  ) %>%
    as_tibble() %>%
    left_join(
      month_tab %>%
        select(
          `..ym..`,
          `..year..`,
          age_band_std,
          `..sex_std..`,
          `..decile_std..`,
          n_m,
          p_m
        ),
      by = c(
        "..ym..",
        "..year..",
        "age_band_std",
        "..sex_std..",
        "..decile_std.."
      )
    )
  
  
  # ---- calculate raw population standardisation weights ----
  # Raw weight for month m and stratum s:
  #   w_raw[m, s] = p_ref[year, s] / p_m[m, s]
  #
  # Strata under-represented in the cohort relative to the external population
  # get weight > 1. Over-represented strata get weight < 1.
  wtab <- month_tab %>%
    left_join(
      pop_tab_join %>%
        select(
          `..year..`,
          age_band_std,
          `..sex_std..`,
          `..decile_std..`,
          pop_n,
          p_ref
        ),
      by = c(
        "..year..",
        "age_band_std",
        "..sex_std..",
        "..decile_std.."
      )
    ) %>%
    mutate(
      p_ref = dplyr::coalesce(p_ref, 0),
      pop_n = dplyr::coalesce(pop_n, 0),
      w_raw = p_ref / p_m
    )
  
  # ---- attach raw weights to individual rows ----
  # Each row receives the raw population weight for its:
  #   month x age band x sex x SIMD quintile.
  out <- df2 %>%
    left_join(
      wtab %>%
        select(
          `..ym..`,
          age_band_std,
          `..sex_std..`,
          `..decile_std..`,
          w_raw
        ),
      by = c(
        "..ym..",
        "age_band_std",
        "..sex_std..",
        "..decile_std.."
      )
    )
  
  # ---- trim and normalise weights ----
  # Step 1: compute a month-specific cap:
  #   cap[m] = weight_cap x median raw weight in month m.
  #
  # Step 2: trim raw weights at that cap:
  #   w_trim0 = min(w_raw, cap[m]).
  #
  # Step 3: normalise the trimmed weights so their mean is 1 within each month:
  #   w_trim = w_trim0 / mean(w_trim0 in month m).
  out <- out %>%
    group_by(`..ym..`) %>%
    mutate(
      `..cap..` = weight_cap * median(w_raw, na.rm = TRUE),
      `..w_trim..` = ifelse(
        is.na(w_raw),
        NA_real_,
        pmin(w_raw, `..cap..`)
      ),
      `..w_trim..` = ifelse(
        is.na(`..w_trim..`),
        NA_real_,
        `..w_trim..` / mean(`..w_trim..`, na.rm = TRUE)
      )
    ) %>%
    ungroup() %>%
    mutate(
      !!weight_col := `..w_trim..`
    )
  
  # ---- quality-assurance summary ----
  # For each month, report:
  #   n     = number of rows with non-missing final weight
  #   max_w = maximum final weight
  #   ess   = effective sample size based on final weights
  qa <- out %>%
    filter(!is.na(.data[[weight_col]])) %>%
    group_by(`..ym..`) %>%
    summarise(
      n = n(),
      max_w = max(.data[[weight_col]], na.rm = TRUE),
      ess = (sum(.data[[weight_col]], na.rm = TRUE)^2) /
        sum(.data[[weight_col]]^2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(`..ym..`)
  
  # ---- weight lookup table ----
  # This is useful for checking the underlying month x stratum proportions
  # and the raw weights before trimming/normalisation.
  weights_table <- wtab %>%
    rename(
      year = `..year..`,
      sex_std = `..sex_std..`,
      decile_std = `..decile_std..`
    ) %>%
    select(
      `..ym..`,
      year,
      age_band_std,
      sex_std,
      decile_std,
      n_m,
      p_m,
      pop_n,
      p_ref,
      w_raw
    ) %>%
    arrange(
      `..ym..`,
      age_band_std,
      sex_std,
      decile_std
    )
  
  # ---- remove temporary calculation columns ----
  out <- out %>%
    select(
      -`..ym..`,
      -`..year..`,
      -`..age_num..`,
      -`..sex_std..`,
      -`..decile_std..`,
      -w_raw,
      -`..cap..`,
      -`..w_trim..`
    )
  
  # ---- return cleaned data and diagnostics ----
  list(
    data = out,
    weights_table = weights_table,
    qa = qa
  )
}


# =============================================================================
# 4) Convenience wrapper: overall + cohort-specific age weights
# =============================================================================
add_all_age_weights_overall_and_subcohort <- function(
    df,
    time_col   = "time",
    age_col    = "age",
    ref_month  = "2012-01",
    age_breaks = AGE_BREAKS_STD,
    age_labels = AGE_LABELS_STD,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = 3.5
) {
  suppressPackageStartupMessages({ library(dplyr) })
  
  # Overall age weights
  res_overall <- add_age_standardisation_weights(
    df = df,
    time_col = time_col,
    age_col = age_col,
    ref_month = ref_month,
    subset_filter = NULL,
    prefix = "w_age",
    age_breaks = age_breaks,
    age_labels = age_labels,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = weight_cap
  )
  df2 <- res_overall$data
  
  # Cohort-specific age weights
  res_fe <- add_age_standardisation_weights(
    df = df2,
    time_col = time_col,
    age_col = age_col,
    ref_month = ref_month,
    subset_filter = quote(dplyr::coalesce(subcohort_FE, 0) == 1),
    prefix = "w_age_FE",
    age_breaks = age_breaks,
    age_labels = age_labels,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = weight_cap
  )
  df2 <- res_fe$data
  
  res_ltc <- add_age_standardisation_weights(
    df = df2,
    time_col = time_col,
    age_col = age_col,
    ref_month = ref_month,
    subset_filter = quote(dplyr::coalesce(subcohort_LTC, 0) == 1),
    prefix = "w_age_LTC",
    age_breaks = age_breaks,
    age_labels = age_labels,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = weight_cap
  )
  df2 <- res_ltc$data
  
  res_yed <- add_age_standardisation_weights(
    df = df2,
    time_col = time_col,
    age_col = age_col,
    ref_month = ref_month,
    subset_filter = quote(dplyr::coalesce(subcohort_YED, 0) == 1),
    prefix = "w_age_YED",
    age_breaks = age_breaks,
    age_labels = age_labels,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = weight_cap
  )
  df2 <- res_yed$data
  
  list(
    data = df2,
    qa_overall = res_overall$qa,
    qa_by_cohort = list(FE = res_fe$qa, LTC = res_ltc$qa, YED = res_yed$qa),
    weights_tables = list(
      Overall = res_overall$weights_table,
      FE = res_fe$weights_table,
      LTC = res_ltc$weights_table,
      YED = res_yed$weights_table
    )
  )
}


# =============================================================================
# 5) Convenience wrapper: add all weights
#    - age-to-reference-month weights
#    - population weights
# =============================================================================
add_all_weights <- function(
    df,
    time_col   = "time",
    age_col    = "age",
    sex_col    = "sexM",
    decile_col = "decile",
    ref_month  = "2012-01",
    pop_csv    = POP_STANDARDISATION_CSV,
    age_breaks = AGE_BREAKS_STD,
    age_labels = AGE_LABELS_STD,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = 3.5
) {
  res_age <- add_all_age_weights_overall_and_subcohort(
    df = df,
    time_col = time_col,
    age_col = age_col,
    ref_month = ref_month,
    age_breaks = age_breaks,
    age_labels = age_labels,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = weight_cap
  )
  
  df2 <- res_age$data
  
  res_pop <- add_population_standardisation_weights(
    df = df2,
    pop_csv = pop_csv,
    time_col = time_col,
    age_col = age_col,
    sex_col = sex_col,
    decile_col = decile_col,
    prefix = "w_pop",
    age_breaks = age_breaks,
    age_labels = age_labels,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = weight_cap
  )
  
  list(
    data = res_pop$data,
    qa_overall = res_age$qa_overall,
    qa_by_cohort = res_age$qa_by_cohort,
    qa_population = res_pop$qa,
    weights_tables = c(
      res_age$weights_tables,
      list(Population = res_pop$weights_table)
    )
  )
}

# =============================================================================
# 6) Dataset builders
# =============================================================================

make_sparra_v3_datasets <- function(
    df_in,
    pop_csv = POP_STANDARDISATION_CSV,
    min_age_included = MIN_AGE_INCLUDED
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyverse)
  })
  
  # --- 1) Keep only variables we need ----------
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
    "age","decile","parkinsons_indicated","MS_indicated","epilepsy_indicated",
    "dementia_indicated", "sexM"
  )
  
  df_v3 <- df_in[, v3_vars[v3_vars %in% names(df_in)], drop = FALSE]
  
  # --- EXCLUSION CRITERIA ----------
  df_v3 <- df_v3[complete.cases(df_v3$target), ]
  df_v3$decile[df_v3$decile == 11] <- NA
  df_v3 <- df_v3[complete.cases(df_v3$decile), ]
  df_v3 <- df_v3[complete.cases(df_v3$age), ]
  df_v3 <- df_v3[df_v3$age >= min_age_included, ]
  # --------------------------------------------
  
  
  # --- add all weights BEFORE cohort-specific feature cuts ---
  df_v3 <- add_all_weights(
    df = df_v3,
    time_col = "time",
    age_col  = "age",
    sex_col  = "sexM",
    decile_col = "decile",
    ref_month = "2012-01",
    pop_csv = pop_csv,
    age_breaks = AGE_BREAKS_STD,
    age_labels = AGE_LABELS_STD,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = 3.5
  )$data
  
  coalesce_counts <- function(.data, cols) {
    mutate(.data, across(any_of(cols), ~ dplyr::coalesce(., 0)))
  }
  
  # --- 2) FE ---------------------------------------------------------------
  fe_counts <- c(
    "num_emergency_admissions","emergency_bed_days","num_alcohol_admissions",
    "num_ae2_attendances","num_elective_admissions",
    "num_outpatient_appointment_general","num_bnf_sections",
    "pis_incontinence","pis_respiratory","pis_cns","pis_infections",
    "pis_endocrine","parkinsons_indicated","numLTCs_resulting_in_admin"
  )
  
  df_FE <- df_v3 %>%
    filter(dplyr::coalesce(subcohort_FE, 0) == 1) %>%
    coalesce_counts(fe_counts) %>%
    mutate(
      age_weight_raw = age,
      decile_weight_raw = decile,
      num_emergency_admissions = cut(num_emergency_admissions,
                                     breaks = c(-Inf,0,1,2,3,4,5,Inf),
                                     labels = c("0","1","2","3","4","5","6+"),
                                     right = TRUE, ordered_result = TRUE),
      emergency_bed_days = cut(emergency_bed_days,
                               breaks = c(-Inf,1,7,14,Inf),
                               labels = c("0-1","2-7","8-14","15+"),
                               right = TRUE, ordered_result = TRUE),
      num_alcohol_admissions = factor(ifelse(num_alcohol_admissions >= 1,"1+","0"),
                                      levels = c("0","1+"), ordered = TRUE),
      num_ae2_attendances = cut(num_ae2_attendances,
                                breaks = c(-Inf,0,3,6,9,Inf),
                                labels = c("0","1-3","4-6","7-9","10+"),
                                right = TRUE, ordered_result = TRUE),
      num_elective_admissions = cut(num_elective_admissions,
                                    breaks = c(-Inf,2,4,6,Inf),
                                    labels = c("0-2","3-4","5-6","7+"),
                                    right = TRUE, ordered_result = TRUE),
      num_outpatient_appointment_general = cut(num_outpatient_appointment_general,
                                               breaks = c(-Inf,0,2,4,Inf),
                                               labels = c("0","1-2","3-4","5+"),
                                               right = TRUE, ordered_result = TRUE),
      num_bnf_sections = cut(num_bnf_sections,
                             breaks = c(-Inf,1,4,7,10,13,16,19,Inf),
                             labels = c("0-1","2-4","5-7","8-10","11-13","14-16","17-19","20+"),
                             right = TRUE, ordered_result = TRUE),
      pis_incontinence = factor(ifelse(pis_incontinence >= 1,"1+","0"),
                                levels = c("0","1+"), ordered = TRUE),
      pis_respiratory  = factor(ifelse(pis_respiratory  >= 8,"8+","0-7"),
                                levels = c("0-7","8+"), ordered = TRUE),
      pis_cns = factor(ifelse(pis_cns >= 8,"8+",
                              ifelse(pis_cns >= 1,"1-7","0")),
                       levels = c("0","1-7","8+"), ordered = TRUE),
      pis_infections = factor(ifelse(pis_infections >= 8,"8+",
                                     ifelse(pis_infections >= 1,"1-7","0")),
                              levels = c("0","1-7","8+"), ordered = TRUE),
      pis_endocrine = factor(ifelse(pis_endocrine >= 8,"8+",
                                    ifelse(pis_endocrine >= 1,"1-7","0")),
                             levels = c("0","1-7","8+"), ordered = TRUE),
      parkinsons_indicated = factor(ifelse(parkinsons_indicated >= 1,"1+","0"),
                                    levels = c("0","1+"), ordered = TRUE),
      numLTCs_resulting_in_admin = cut(numLTCs_resulting_in_admin,
                                       breaks = c(-Inf,0,1,2,3,4,5,Inf),
                                       labels = c("0","1","2","3","4","5","6+"),
                                       right = TRUE, ordered_result = TRUE),
      age = cut(age,
                breaks = c(74,79,84,89,Inf),
                labels = c("75-79","80-84","85-89","90+"),
                right = TRUE, ordered_result = TRUE),
      decile = cut(decile,
                   breaks = c(-Inf,2,4,6,8,10),
                   labels = c("1","2","3","4","5"),
                   right = TRUE, ordered_result = TRUE),
      int_age__num_emergency_admissions =
        interaction(age, num_emergency_admissions, sep="__", drop=TRUE),
      int_age__numLTCs_resulting_in_admin =
        interaction(age, numLTCs_resulting_in_admin, sep="__", drop=TRUE),
      int_num_bnf_sections__pis_respiratory =
        interaction(num_bnf_sections, pis_respiratory, sep="__", drop=TRUE)
    ) %>%
    select(any_of(c(
      "id","time","target",
      "num_emergency_admissions","emergency_bed_days","num_alcohol_admissions",
      "num_ae2_attendances","num_elective_admissions",
      "num_outpatient_appointment_general","num_bnf_sections",
      "pis_incontinence","pis_respiratory","pis_cns","pis_infections","pis_endocrine",
      "parkinsons_indicated","numLTCs_resulting_in_admin",
      "age","decile","age_weight_raw","decile_weight_raw",
      "int_age__num_emergency_admissions","int_age__numLTCs_resulting_in_admin",
      "int_num_bnf_sections__pis_respiratory", "sexM",
      "age_band_std",
      "w_age_raw","w_age_stab","w_age_trim",
      "w_age_FE_raw","w_age_FE_stab","w_age_FE_trim",
      "w_pop_raw","w_pop_stab","w_pop_trim"
    )))
  
  # --- 3) LTC --------------------------------------------------------------
  ltc_counts <- c(
    "num_emergency_admissions","emergency_bed_days",
    "num_alcohol_drug_admissions","num_ae2_attendances",
    "num_daycase_admissions","num_elective_admissions","elective_bed_days",
    "num_psych_admissions","num_outpatient_appointment_general","num_bnf_sections",
    "pis_infections","pis_sub_depend","pis_dementia","pis_corticosteroids",
    "pis_fluids","pis_nutrition","pis_vitamins","pis_bandages","pis_catheters","pis_stoma",
    "numLTCs_resulting_in_admin","epilepsy_indicated","MS_indicated","parkinsons_indicated",
    "decile"
  )
  
  df_LTC <- df_v3 %>%
    filter(dplyr::coalesce(subcohort_LTC, 0) == 1) %>%
    coalesce_counts(ltc_counts) %>%
    mutate(
      age_weight_raw = age,
      decile_weight_raw = decile,
      num_emergency_admissions = cut(num_emergency_admissions,
                                     breaks = c(-Inf,0,1,2,3,4,5,Inf),
                                     labels = c("0","1","2","3","4","5","6+"),
                                     right = TRUE, ordered_result = TRUE),
      emergency_bed_days = cut(emergency_bed_days,
                               breaks = c(-Inf,1,7,14,21,Inf),
                               labels = c("0-1","2-7","8-14","15-21","22+"),
                               right = TRUE, ordered_result = TRUE),
      num_alcohol_drug_admissions = cut(num_alcohol_drug_admissions,
                                        breaks = c(-Inf,0,1,2,Inf),
                                        labels = c("0","1","2","3+"),
                                        right = TRUE, ordered_result = TRUE),
      num_ae2_attendances = cut(num_ae2_attendances,
                                breaks = c(-Inf,0,3,6,9,Inf),
                                labels = c("0","1-3","4-6","7-9","10+"),
                                right = TRUE, ordered_result = TRUE),
      num_daycase_admissions = cut(num_daycase_admissions,
                                   breaks = c(-Inf,0,1,2,Inf),
                                   labels = c("0","1","2","3+"),
                                   right = TRUE, ordered_result = TRUE),
      num_elective_admissions = cut(num_elective_admissions,
                                    breaks = c(-Inf,0,1,2,Inf),
                                    labels = c("0","1","2","3+"),
                                    right = TRUE, ordered_result = TRUE),
      elective_bed_days = cut(elective_bed_days,
                              breaks = c(-Inf,28,Inf),
                              labels = c("0-28","29+"),
                              right = TRUE, ordered_result = TRUE),
      num_psych_admissions = factor(ifelse(num_psych_admissions >= 1, "1+", "0"),
                                    levels = c("0","1+"), ordered = TRUE),
      num_outpatient_appointment_general = cut(num_outpatient_appointment_general,
                                               breaks = c(-Inf,0,2,4,6,Inf),
                                               labels = c("0","1-2","3-4","5-6","7+"),
                                               right = TRUE, ordered_result = TRUE),
      numLTCs_resulting_in_admin = cut(numLTCs_resulting_in_admin,
                                       breaks = c(-Inf,0,1,2,3,4,5,Inf),
                                       labels = c("0","1","2","3","4","5","6+"),
                                       right = TRUE, ordered_result = TRUE),
      num_bnf_sections = cut(num_bnf_sections,
                             breaks = c(-Inf,1,4,7,10,13,Inf),
                             labels = c("0-1","2-4","5-7","8-10","11-13","14+"),
                             right = TRUE, ordered_result = TRUE),
      epilepsy_indicated = factor(dplyr::case_when(
        epilepsy_indicated <= 0 ~ "0",
        epilepsy_indicated == 1 ~ "1",
        epilepsy_indicated >= 2 ~ "2"
      ), levels = c("0","1","2"), ordered = TRUE),
      MS_indicated = factor(ifelse(MS_indicated >= 2, "2", "0-1"),
                            levels = c("0-1","2"), ordered = TRUE),
      parkinsons_indicated = factor(ifelse(parkinsons_indicated >= 2, "2", "0-1"),
                                    levels = c("0-1","2"), ordered = TRUE),
      pis_infections = factor(ifelse(pis_infections >= 8, "8+", "0-7"),
                              levels = c("0-7","8+"), ordered = TRUE),
      pis_sub_depend = factor(ifelse(pis_sub_depend >= 1, "1+", "0"),
                              levels = c("0","1+"), ordered = TRUE),
      pis_dementia = factor(ifelse(pis_dementia >= 1, "1+", "0"),
                            levels = c("0","1+"), ordered = TRUE),
      pis_corticosteroids = factor(ifelse(pis_corticosteroids >= 1, "1+", "0"),
                                   levels = c("0","1+"), ordered = TRUE),
      pis_fluids = factor(ifelse(pis_fluids >= 1, "1+", "0"),
                          levels = c("0","1+"), ordered = TRUE),
      pis_nutrition = factor(ifelse(pis_nutrition >= 1, "1+", "0"),
                             levels = c("0","1+"), ordered = TRUE),
      pis_vitamins = factor(ifelse(pis_vitamins >= 2, "2", "0-1"),
                            levels = c("0-1","2"), ordered = TRUE),
      pis_bandages = factor(ifelse(pis_bandages >= 2, "2", "0-1"),
                            levels = c("0-1","2"), ordered = TRUE),
      pis_catheters = factor(ifelse(pis_catheters >= 1, "1+", "0"),
                             levels = c("0","1+"), ordered = TRUE),
      pis_stoma = factor(ifelse(pis_stoma >= 1, "1+", "0"),
                         levels = c("0","1+"), ordered = TRUE),
      decile = cut(decile,
                   breaks = c(-Inf,2,4,6,8,10),
                   labels = c("1","2","3","4","5"),
                   right = TRUE, ordered_result = TRUE),
      int_num_alcohol_drug_admissions__pis_sub_depend =
        interaction(num_alcohol_drug_admissions, pis_sub_depend, sep="__", drop=TRUE),
      int_num_emergency_admissions__numLTCs_resulting_in_admin =
        interaction(num_emergency_admissions, numLTCs_resulting_in_admin, sep="__", drop=TRUE),
      int_num_psych_admissions__pis_dementia =
        interaction(num_psych_admissions, pis_dementia, sep="__", drop=TRUE)
    ) %>%
    select(any_of(c(
      "id","time","target",
      "num_emergency_admissions","emergency_bed_days",
      "num_alcohol_drug_admissions","num_ae2_attendances",
      "num_daycase_admissions","num_elective_admissions","elective_bed_days",
      "num_psych_admissions","num_outpatient_appointment_general",
      "numLTCs_resulting_in_admin","num_bnf_sections",
      "epilepsy_indicated","MS_indicated","parkinsons_indicated",
      "pis_infections","pis_sub_depend","pis_dementia","pis_corticosteroids",
      "pis_fluids","pis_nutrition","pis_vitamins","pis_bandages","pis_catheters","pis_stoma",
      "decile","age_weight_raw","decile_weight_raw",
      "int_num_alcohol_drug_admissions__pis_sub_depend",
      "int_num_emergency_admissions__numLTCs_resulting_in_admin",
      "int_num_psych_admissions__pis_dementia", "sexM",
      "age_band_std",
      "w_age_raw","w_age_stab","w_age_trim",
      "w_age_LTC_raw","w_age_LTC_stab","w_age_LTC_trim",
      "w_pop_raw","w_pop_stab","w_pop_trim"
    )))
  
  # --- 4) YED --------------------------------------------------------------
  yed_counts <- c(
    "num_alcohol_drug_admissions","num_emergency_admissions","emergency_bed_days",
    "num_ae2_attendances","num_psych_admissions","num_elective_admissions",
    "num_outpatient_appointment_general","num_outpatient_appointment_psych",
    "numLTCs_resulting_in_admin","num_bnf_sections",
    "pis_gut_motility","pis_antisecretory","pis_intestinal",
    "pis_antifibrinolytic","pis_anticoagulant","pis_stoma","pis_mucolytics",
    "pis_diabetes","pis_corticosteroids","pis_fluids","pis_vitamins",
    "pis_cns","dementia_indicated"
  )
  
  df_YED <- df_v3 %>%
    filter(dplyr::coalesce(subcohort_YED, 0) == 1) %>%
    coalesce_counts(yed_counts) %>%
    mutate(
      age_weight_raw = age,
      decile_weight_raw = decile,
      num_alcohol_drug_admissions = cut(num_alcohol_drug_admissions,
                                        breaks = c(-Inf, 0, 1, Inf),
                                        labels = c("0","1","2+"),
                                        right = TRUE, ordered_result = TRUE),
      num_emergency_admissions = cut(num_emergency_admissions,
                                     breaks = c(-Inf,0,1,2,3,4,5,Inf),
                                     labels = c("0","1","2","3","4","5","6+"),
                                     right = TRUE, ordered_result = TRUE),
      emergency_bed_days = cut(emergency_bed_days,
                               breaks = c(-Inf,1,7,14,21,Inf),
                               labels = c("0-1","2-7","8-14","15-21","22+"),
                               right = TRUE, ordered_result = TRUE),
      num_ae2_attendances = cut(num_ae2_attendances,
                                breaks = c(-Inf,3,6,9,Inf),
                                labels = c("0-3","4-6","7-9","10+"),
                                right = TRUE, ordered_result = TRUE),
      num_psych_admissions = factor(ifelse(num_psych_admissions >= 1, "1+", "0"),
                                    levels = c("0","1+"), ordered = TRUE),
      num_elective_admissions = cut(num_elective_admissions,
                                    breaks = c(-Inf, 2, 4, 6, Inf),
                                    labels = c("0-2","3-4","5-6","7+"),
                                    right = TRUE, ordered_result = TRUE),
      num_outpatient_appointment_general = cut(num_outpatient_appointment_general,
                                               breaks = c(-Inf,0,2,4,Inf),
                                               labels = c("0","1-2","3-4","5+"),
                                               right = TRUE, ordered_result = TRUE),
      num_outpatient_appointment_psych =
        factor(ifelse(num_outpatient_appointment_psych >= 1, "1+", "0"),
               levels = c("0","1+"), ordered = TRUE),
      numLTCs_resulting_in_admin = cut(numLTCs_resulting_in_admin,
                                       breaks = c(-Inf,0,1,2,Inf),
                                       labels = c("0","1","2","3+"),
                                       right = TRUE, ordered_result = TRUE),
      num_bnf_sections = cut(num_bnf_sections,
                             breaks = c(-Inf,1,4,7,10,13,Inf),
                             labels = c("0-1","2-4","5-7","8-10","11-13","14+"),
                             right = TRUE, ordered_result = TRUE),
      pis_gut_motility     = factor(ifelse(pis_gut_motility     >= 1, "1+", "0"), levels = c("0","1+"), ordered = TRUE),
      pis_antisecretory    = factor(ifelse(pis_antisecretory    >= 1, "1+", "0"), levels = c("0","1+"), ordered = TRUE),
      pis_intestinal       = factor(ifelse(pis_intestinal       >= 1, "1+", "0"), levels = c("0","1+"), ordered = TRUE),
      pis_antifibrinolytic = factor(ifelse(pis_antifibrinolytic >= 1, "1+", "0"), levels = c("0","1+"), ordered = TRUE),
      pis_anticoagulant    = factor(ifelse(pis_anticoagulant    >= 1, "1+", "0"), levels = c("0","1+"), ordered = TRUE),
      pis_stoma            = factor(ifelse(pis_stoma            >= 1, "1+", "0"), levels = c("0","1+"), ordered = TRUE),
      pis_mucolytics       = factor(ifelse(pis_mucolytics       >= 1, "1+", "0"), levels = c("0","1+"), ordered = TRUE),
      pis_diabetes         = factor(ifelse(pis_diabetes         >= 1, "1+", "0"), levels = c("0","1+"), ordered = TRUE),
      pis_corticosteroids  = factor(ifelse(pis_corticosteroids  >= 1, "1+", "0"), levels = c("0","1+"), ordered = TRUE),
      pis_fluids           = factor(ifelse(pis_fluids           >= 1, "1+", "0"), levels = c("0","1+"), ordered = TRUE),
      pis_vitamins         = factor(ifelse(pis_vitamins         >= 1, "1+", "0"), levels = c("0","1+"), ordered = TRUE),
      pis_cns = factor(
        ifelse(pis_cns >= 8, "8+", ifelse(pis_cns >= 1, "1-7", "0")),
        levels = c("0","1-7","8+"), ordered = TRUE
      ),
      dementia_indicated   = factor(ifelse(dementia_indicated   >= 1, "1+", "0"),
                                    levels = c("0","1+"), ordered = TRUE),
      int_num_psych_admissions__num_alcohol_drug_admissions =
        interaction(num_psych_admissions, num_alcohol_drug_admissions, sep="__", drop=TRUE),
      int_numLTCs_resulting_in_admin__num_bnf_sections =
        interaction(numLTCs_resulting_in_admin, num_bnf_sections, sep="__", drop=TRUE),
      int_num_alcohol_drug_admissions__pis_cns =
        interaction(num_alcohol_drug_admissions, pis_cns, sep="__", drop=TRUE),
      int_num_bnf_sections__pis_vitamins =
        interaction(num_bnf_sections, pis_vitamins, sep="__", drop=TRUE)
    ) %>%
    select(any_of(c(
      "id","time","target",
      "num_alcohol_drug_admissions",
      "num_emergency_admissions","emergency_bed_days",
      "num_ae2_attendances","num_psych_admissions",
      "num_elective_admissions",
      "num_outpatient_appointment_general","num_outpatient_appointment_psych",
      "numLTCs_resulting_in_admin","num_bnf_sections",
      "pis_gut_motility","pis_antisecretory","pis_intestinal",
      "pis_antifibrinolytic","pis_anticoagulant","pis_stoma","pis_mucolytics",
      "pis_diabetes","pis_corticosteroids","pis_fluids","pis_vitamins",
      "pis_cns","dementia_indicated","age_weight_raw","decile_weight_raw",
      "int_num_psych_admissions__num_alcohol_drug_admissions",
      "int_numLTCs_resulting_in_admin__num_bnf_sections",
      "int_num_alcohol_drug_admissions__pis_cns",
      "int_num_bnf_sections__pis_vitamins", "sexM",
      "age_band_std",
      "w_age_raw","w_age_stab","w_age_trim",
      "w_age_YED_raw","w_age_YED_stab","w_age_YED_trim",
      "w_pop_raw","w_pop_stab","w_pop_trim"
    )))
  
  overall_vars <- union(
    union(names(df_FE), names(df_LTC)),
    names(df_YED)
  )

  df_overall <- df_v3 %>%
    mutate(
      age_weight_raw = age,
      decile_weight_raw = decile
    ) %>%
    select(any_of(overall_vars))

  list(FE = df_FE, LTC = df_LTC, YED = df_YED, Overall = df_overall)
}


make_sparra_v3_datasets_raw <- function(
    df_in,
    pop_csv = POP_STANDARDISATION_CSV,
    min_age_included = MIN_AGE_INCLUDED
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyverse)
  })
  
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
    "age","decile","parkinsons_indicated","MS_indicated","epilepsy_indicated",
    "dementia_indicated",
    "first_event_reason", "first_event_time",
    "emerg_reason", "emerg_first_time", "sexM"
  )
  
  df_v3 <- df_in[, v3_vars[v3_vars %in% names(df_in)], drop = FALSE]
  
  # --- EXCLUSION CRITERIA ----------
  df_v3 <- df_v3[complete.cases(df_v3$target), ]
  df_v3$decile[df_v3$decile == 11] <- NA
  df_v3 <- df_v3[complete.cases(df_v3$decile), ]
  df_v3 <- df_v3[complete.cases(df_v3$age), ]
  df_v3 <- df_v3[df_v3$age >= min_age_included, ]
  # --------------------------------------------
  
  df_v3 <- add_all_weights(
    df = df_v3,
    time_col = "time",
    age_col  = "age",
    sex_col  = "sexM",
    decile_col = "decile",
    ref_month = "2012-01",
    pop_csv = pop_csv,
    age_breaks = AGE_BREAKS_STD,
    age_labels = AGE_LABELS_STD,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = 3.5
  )$data
  
  df_FE  <- df_v3 %>% filter(dplyr::coalesce(subcohort_FE, 0) == 1)
  df_LTC <- df_v3 %>% filter(dplyr::coalesce(subcohort_LTC, 0) == 1)
  df_YED <- df_v3 %>% filter(dplyr::coalesce(subcohort_YED, 0) == 1)
  
  df_overall <- df_v3

  list(FE = df_FE, LTC = df_LTC, YED = df_YED, Overall = df_overall)
}