# =============================================================================
# standardisationHelpers.R
#
# Author       : Louis Chislett
#
# Purpose :
# Prepare and reuse standardisation inputs for the SPARRA v3 evaluation pipeline,
# and recalculate bootstrap weights for reference-weighted and population-
# weighted model summaries.
#
# This script:
# (1) Prepares a population standardisation lookup
# - Reads the population CSV once
# - Standardises sex, SIMD decile/quintile, and age-band fields
# - Calculates year-specific population proportions for each age x sex x SIMD
#   cell used in population weighting
#
# (2) Caches the population lookup
# - Avoids rereading and rebuilding the same population table during repeated
#   bootstrap calculations
# - Uses the population file path and age-band settings as the cache key
#
# (3) Recalculates reference-month bootstrap weights
# - Combines the fixed reference month with one bootstrapped evaluation month
# - Calls the existing age-standardisation helper
# - Returns bootstrap-row reference weights aligned to the sampled row positions
#
# (4) Recalculates population bootstrap weights
# - Uses the cached population lookup rather than rereading the CSV
# - Recreates population weights for a bootstrapped evaluation month
# - Returns bootstrap-row population weights aligned to the sampled row positions
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(tidyr)
})


# Cache for prepared population standardisation tables.
#
# The cache is intentionally file/settings-specific, so the same population CSV
# can be reused across many bootstrap replicates without being read repeatedly.
.pop_standardisation_cache <- new.env(parent = emptyenv())


# -----------------------------------------------------------------------------
# Population standardisation lookup
# -----------------------------------------------------------------------------

# Read the population standardisation CSV and convert it into the lookup table
# needed for population weighting.
#
# Output columns include:
#   year
#   age_band_std
#   ..sex_std..
#   ..decile_std..
#   pop_n
#   p_ref
#
# p_ref is the population proportion within a year for each age x sex x SIMD
# cell. This is later compared with the corresponding cohort proportion p_m.
prepare_population_standardisation_lookup <- function(pop_csv,
                                                      age_breaks = AGE_BREAKS_STD,
                                                      age_labels = AGE_LABELS_STD,
                                                      min_age_included = MIN_AGE_INCLUDED) {
  needed_pop_cols <- c("year", "sex", "simd_decile", "age", "population")
  
  pop_raw <- readr::read_csv(pop_csv, show_col_types = FALSE)
  
  missing_pop_cols <- setdiff(needed_pop_cols, names(pop_raw))
  if (length(missing_pop_cols) > 0) {
    stop(sprintf(
      "Population CSV is missing required column(s): %s",
      paste(missing_pop_cols, collapse = ", ")
    ))
  }
  
  pop_raw %>%
    mutate(
      year = as.integer(year),
      sex_std = standardise_sex_from_population(sex),
      simd_decile = as.integer(simd_decile),
      simd_quintile = decile_to_quintile(simd_decile),
      age_num = parse_population_age_to_numeric(age),
      population = as.numeric(population)
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
    summarise(pop_n = sum(population, na.rm = TRUE), .groups = "drop") %>%
    group_by(year) %>%
    mutate(p_ref = pop_n / sum(pop_n)) %>%
    ungroup() %>%
    rename(
      `..sex_std..` = sex_std,
      `..decile_std..` = simd_quintile
    )
}


# Return a cached population standardisation lookup where possible.
# This is useful because bootstrap code can request the same population table many
# times. Rebuilding it for every replicate would be unnecessary and slow.
get_population_standardisation_lookup_cached <- function(pop_csv,
                                                         age_breaks = AGE_BREAKS_STD,
                                                         age_labels = AGE_LABELS_STD,
                                                         min_age_included = MIN_AGE_INCLUDED) {
  key <- paste(
    normalizePath(pop_csv, winslash = "/", mustWork = FALSE),
    paste(age_breaks, collapse = ","),
    paste(age_labels, collapse = ","),
    min_age_included,
    sep = "||"
  )
  
  if (!exists(key, envir = .pop_standardisation_cache, inherits = FALSE)) {
    .pop_standardisation_cache[[key]] <- prepare_population_standardisation_lookup(
      pop_csv = pop_csv,
      age_breaks = age_breaks,
      age_labels = age_labels,
      min_age_included = min_age_included
    )
  }
  
  .pop_standardisation_cache[[key]]
}


# -----------------------------------------------------------------------------
# Bootstrap reference weights
# -----------------------------------------------------------------------------

# Recalculate reference-month age-standardised weights for one bootstrap sample.
#
# month_df is the evaluation-month data.
# ref_df is the fixed reference-month data, usually the train/reference month.
# idx is the sampled row index vector for the bootstrap replicate.
#
# The function returns one row per usable sampled bootstrap row:
#   idx      = original row index in month_df
#   boot_pos = position in the bootstrap sample
#   w_ref    = recalculated reference weight
#
# NULL is returned for empty or unusable bootstrap samples. This is kept because
# upstream bootstrap code can skip invalid replicates cleanly.
.recalc_reference_weights_for_idx <- function(month_df,
                                              ref_df,
                                              idx,
                                              ref_month = "2012-01",
                                              age_col = "age_weight_raw",
                                              ref_prefix = "boot_w_ref",
                                              weight_cap = 3.5) {
  if (!length(idx) || !nrow(month_df) || !nrow(ref_df)) return(NULL)
  if (!all(c("time", "target", age_col) %in% names(month_df))) return(NULL)
  if (!all(c("time", age_col) %in% names(ref_df))) return(NULL)
  
  # Keep the sampled row position and original row index so weights can be
  # returned in bootstrap-sample order after recalculation.
  month_boot <- month_df[idx, , drop = FALSE]
  month_boot$..boot_pos.. <- seq_len(nrow(month_boot))
  month_boot$..orig_idx.. <- as.integer(idx)
  
  # The reference month is fixed, not resampled. It is included only so the
  # standardisation helper can calculate reference-based weights.
  ref_fixed <- ref_df
  ref_fixed$..boot_pos.. <- NA_integer_
  ref_fixed$..orig_idx.. <- NA_integer_
  
  # Align columns before binding the reference rows and bootstrapped month rows.
  keep_cols <- union(names(ref_fixed), names(month_boot))
  ref_fixed <- ref_fixed[, intersect(keep_cols, names(ref_fixed)), drop = FALSE]
  month_boot <- month_boot[, intersect(keep_cols, names(month_boot)), drop = FALSE]
  
  boot_dat <- bind_rows(ref_fixed, month_boot)
  
  weighted_boot <- add_age_standardisation_weights(
    df = boot_dat,
    time_col = "time",
    age_col = age_col,
    ref_month = ref_month,
    subset_filter = NULL,
    prefix = ref_prefix,
    age_breaks = AGE_BREAKS_STD,
    age_labels = AGE_LABELS_STD,
    min_age_included = MIN_AGE_INCLUDED,
    weight_cap = weight_cap
  )$data
  
  weight_col <- paste0(ref_prefix, "_trim")
  if (!weight_col %in% names(weighted_boot)) return(NULL)
  
  out <- weighted_boot %>%
    filter(!is.na(.data[["..boot_pos.."]]), !is.na(.data[["..orig_idx.."]])) %>%
    arrange(.data[["..boot_pos.."]]) %>%
    transmute(
      idx = as.integer(.data[["..orig_idx.."]]),
      boot_pos = as.integer(.data[["..boot_pos.."]]),
      w_ref = sanitize_weights(.data[[weight_col]])
    ) %>%
    filter(is.finite(.data[["w_ref"]]), !is.na(.data[["w_ref"]]), .data[["w_ref"]] > 0)
  
  if (!nrow(out)) return(NULL)
  
  out
}


# -----------------------------------------------------------------------------
# Bootstrap population weights
# -----------------------------------------------------------------------------

# Recalculate population-standardised weights for one bootstrap sample.
#
# This is the cached/bootstrap equivalent of add_population_standardisation_weights().
# It avoids reading the population CSV inside every bootstrap replicate by using
# the prepared pop_lookup table from get_population_standardisation_lookup_cached().
#
# The calculation mirrors the main cohort-builder population weighting:
#   - classify sampled rows by month, year, age band, sex, and SIMD quintile,
#   - calculate cohort cell proportions p_m,
#   - join the population proportions p_ref,
#   - calculate raw weights p_ref / p_m,
#   - trim at weight_cap * median(raw weight),
#   - normalise weights to mean 1 within month.
#
# The function returns one row per usable sampled bootstrap row:
#   idx      = original row index in month_df
#   boot_pos = position in the bootstrap sample
#   w_pop    = recalculated population weight
.recalc_population_weights_for_idx <- function(month_df,
                                               pop_lookup,
                                               idx,
                                               age_col = "age_weight_raw",
                                               sex_col = "sexM",
                                               decile_col = "decile_weight_raw",
                                               weight_cap = 3.5) {
  if (!length(idx) || !nrow(month_df) || is.null(pop_lookup) || !nrow(pop_lookup)) return(NULL)
  if (!all(c("time", "target", age_col, sex_col, decile_col) %in% names(month_df))) return(NULL)
  
  # Keep sampled row position and original row index for later alignment.
  month_boot <- month_df[idx, , drop = FALSE]
  month_boot$..boot_pos.. <- seq_len(nrow(month_boot))
  month_boot$..orig_idx.. <- as.integer(idx)
  
  # Standardise cohort-side fields so they can join to the population lookup.
  cohort_boot <- month_boot %>%
    mutate(
      `..ym..` = to_ym(.data[["time"]]),
      `..year..` = as.integer(substr(`..ym..`, 1, 4)),
      `..age_num..` = as.numeric(.data[[age_col]]),
      `..sex_std..` = standardise_sex_from_cohort(.data[[sex_col]]),
      `..decile_std..` = decile_to_quintile(.data[[decile_col]]),
      age_band_std = cut(
        `..age_num..`,
        breaks = AGE_BREAKS_STD,
        labels = AGE_LABELS_STD,
        right = FALSE,
        ordered_result = TRUE
      )
    ) %>%
    filter(
      !is.na(`..ym..`),
      !is.na(`..year..`),
      !is.na(`..age_num..`),
      `..age_num..` >= MIN_AGE_INCLUDED,
      !is.na(age_band_std),
      !is.na(`..sex_std..`),
      !is.na(`..decile_std..`)
    ) %>%
    mutate(
      age_band_std = collapse_age_band_for_population_weights(age_band_std)
    )
  
  if (!nrow(cohort_boot)) return(NULL)
  
  # Check that the population lookup contains every year present in the bootstrap
  # sample. If not, this replicate cannot receive valid population weights.
  month_years <- cohort_boot %>% distinct(`..ym..`, `..year..`)
  missing_years <- setdiff(
    sort(unique(month_years$`..year..`)),
    sort(unique(pop_lookup$year))
  )
  
  if (length(missing_years) > 0) return(NULL)
  
  # Cohort-side cell counts and proportions for the bootstrapped month.
  month_tab <- cohort_boot %>%
    count(`..ym..`, `..year..`, age_band_std, `..sex_std..`, `..decile_std..`, name = "n_m") %>%
    group_by(`..ym..`) %>%
    mutate(p_m = n_m / sum(n_m)) %>%
    ungroup()
  
  # Build the full set of population cells expected for each month/year, then
  # join the cohort proportions onto it. This identifies missing cohort cells
  # where the population has non-zero mass.
  full_grid <- purrr::map_dfr(sort(unique(month_years$`..year..`)), function(.yr) {
    my <- month_years %>% filter(`..year..` == .yr)
    py <- pop_lookup %>% filter(year == .yr)
    
    if (!nrow(my) || !nrow(py)) return(tibble())
    
    tidyr::crossing(my, py %>% select(-year))
  }) %>%
    left_join(
      month_tab %>%
        select(`..ym..`, age_band_std, `..sex_std..`, `..decile_std..`, n_m, p_m),
      by = c("..ym..", "age_band_std", "..sex_std..", "..decile_std..")
    )
  
  if (!nrow(full_grid)) return(NULL)
  
  # A replicate is unusable if a population cell has positive population mass but
  # no sampled cohort rows. That would require an infinite population weight.
  has_missing_positive_population_cell <- any(
    full_grid$p_ref > 0 &
      (is.na(full_grid$n_m) | is.na(full_grid$p_m) | full_grid$p_m == 0)
  )
  
  if (has_missing_positive_population_cell) return(NULL)
  
  # Calculate raw population weights for the cells observed in the bootstrapped
  # cohort sample.
  wtab <- month_tab %>%
    left_join(
      pop_lookup %>%
        select(year, age_band_std, `..sex_std..`, `..decile_std..`, pop_n, p_ref),
      by = c("..year.." = "year", "age_band_std", "..sex_std..", "..decile_std..")
    ) %>%
    mutate(
      p_ref = dplyr::coalesce(p_ref, 0),
      pop_n = dplyr::coalesce(pop_n, 0),
      w_raw = p_ref / p_m
    )
  
  # Join cell weights back onto sampled rows, trim extreme weights, normalise to
  # mean 1 within month, and return the weights in bootstrap-sample order.
  weighted_boot <- cohort_boot %>%
    left_join(
      wtab %>% select(`..ym..`, age_band_std, `..sex_std..`, `..decile_std..`, w_raw),
      by = c("..ym..", "age_band_std", "..sex_std..", "..decile_std..")
    ) %>%
    group_by(`..ym..`) %>%
    mutate(
      `..cap..` = weight_cap * median(w_raw, na.rm = TRUE),
      w_pop = ifelse(
        is.na(w_raw),
        NA_real_,
        pmin(w_raw, `..cap..`)
      ),
      w_pop = ifelse(
        is.na(w_pop),
        NA_real_,
        w_pop / mean(w_pop, na.rm = TRUE)
      )
    ) %>%
    ungroup() %>%
    transmute(
      idx = as.integer(.data[["..orig_idx.."]]),
      boot_pos = as.integer(.data[["..boot_pos.."]]),
      w_pop = sanitize_weights(.data[["w_pop"]])
    ) %>%
    filter(is.finite(.data[["w_pop"]]), !is.na(.data[["w_pop"]]), .data[["w_pop"]] > 0) %>%
    arrange(.data[["boot_pos"]])
  
  if (!nrow(weighted_boot)) return(NULL)
  
  weighted_boot
}
