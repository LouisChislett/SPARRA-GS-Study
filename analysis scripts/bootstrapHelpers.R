# =============================================================================
# bootstrapHelpers.R
#
# Author       : Louis Chislett
#
# Purpose :
# Create shared bootstrap replicate objects for the SPARRA v3 fixed odds-ratio
# evaluation pipeline.
#
# This script:
# (1) Draws stratified bootstrap row samples for one evaluation month
# - Samples outcome-positive and outcome-negative rows separately
# - Preserves the original class balance within each bootstrap replicate
# - Returns no sample when a valid two-class bootstrap cannot be drawn
#
# (2) Recalculates reference-month weights for each bootstrap sample
# - Uses the supplied model-specific reference cohort, usually 2012-01
# - Produces reference-weighted bootstrap rows for FE, LTC, YED, and COMBINED
#
# (3) Optionally recalculates population weights for each bootstrap sample
# - Used for the COMBINED population-weighted analysis
# - Skipped for FE, LTC, and YED because we cannot calculate weights for those cohorts
#
# (4) Returns one shared replicate object per successful bootstrap sample
# - The same replicate objects are reused by metrics and SHAP summaries
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})


# -----------------------------------------------------------------------------
# Row sampling
# -----------------------------------------------------------------------------

# Draw one stratified bootstrap sample from a binary outcome vector.
#
# Positive and negative rows are sampled separately, each with replacement and
# each to its original class size. This keeps the bootstrap class balance stable.
#
# Returns integer(0) when the month cannot support stratified resampling, for
# example if there is only one outcome class or fewer than two rows in either
# class.
.stratified_uniform_sample_idx <- function(y) {
  y <- as.integer(y)
  
  if (!length(y) || length(unique(stats::na.omit(y))) < 2L) {
    return(integer(0))
  }
  
  idx_pos <- which(y == 1L)
  idx_neg <- which(y == 0L)
  
  if (length(idx_pos) < 2L || length(idx_neg) < 2L) {
    return(integer(0))
  }
  
  c(
    sample(idx_pos, size = length(idx_pos), replace = TRUE),
    sample(idx_neg, size = length(idx_neg), replace = TRUE)
  )
}


# -----------------------------------------------------------------------------
# Shared bootstrap replicate construction
# -----------------------------------------------------------------------------

# Build shared bootstrap replicate objects for one model/month.
#
# month_df is the evaluation-month data and must contain:
#   target
#   time
#   age_weight_raw
# plus sexM and decile_weight_raw when include_population = TRUE.
#
# ref_df is the fixed reference-month cohort used for reference-weight
# recalculation. For FE/LTC/YED, this should be the corresponding model-specific
# cohort in the reference month. For COMBINED, this should be the Overall
# reference-month population.
#
# pop_lookup is the cached population standardisation lookup. It is only used
# when include_population = TRUE.
#
# Each successful replicate contains:
#   idx                    original row indices sampled from month_df
#   w_ref                  recalculated reference weights
#   w_pop                  recalculated population weights, or NA for FE/LTC/YED
#   bootstrap_rep          requested bootstrap replicate number
#   n_sampled              number of rows sampled before weight filtering
#   n_used                 number of rows retained after weight filtering
#   has_population_weights whether w_pop was calculated
make_shared_bootstrap_reps <- function(month_df,
                                       ref_df,
                                       pop_lookup,
                                       boot_n = 100,
                                       boot_seed = 123,
                                       ref_month = "2012-01",
                                       age_col = "age_weight_raw",
                                       sex_col = "sexM",
                                       decile_col = "decile_weight_raw",
                                       ref_prefix = "boot_w_ref",
                                       weight_cap = 3.5,
                                       include_population = TRUE) {
  if (missing(month_df) || is.null(month_df) || !nrow(month_df)) {
    return(vector("list", 0L))
  }
  
  
  if (is.null(ref_df) || !nrow(ref_df)) {
    return(vector("list", 0L))
  }
  
  boot_n <- as.integer(boot_n)[1]
  
  
  y <- as.integer(month_df$target)
  
  # check that y has two valid values - if it does not then a bootstrapping is skipped
  if (length(unique(stats::na.omit(y))) < 2L) {
    return(vector("list", 0L))
  }
  
  # Set a deterministic seed for this model/month, then restore the caller's RNG
  # state when the function exits so this helper does not affect unrelated code.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  set.seed(as.integer(boot_seed)[1])
  
  reps <- vector("list", boot_n)
  used <- 0L
  
  for (b in seq_len(boot_n)) {
    idx0 <- .stratified_uniform_sample_idx(y)
    
    if (!length(idx0)) {
      next
    }
    
    # Recalculate reference weights for the sampled rows using the supplied
    # reference-month cohort. Invalid samples are skipped.
    ref_b <- .recalc_reference_weights_for_idx(
      month_df = month_df,
      ref_df = ref_df,
      idx = idx0,
      ref_month = ref_month,
      age_col = age_col,
      ref_prefix = paste0(ref_prefix, "_", b),
      weight_cap = weight_cap
    )
    
    if (is.null(ref_b) || !nrow(ref_b)) {
      next
    }
    
    # Population weights are only needed for the Overall/COMBINED population-
    # weighted analysis. FE/LTC/YED are restricted subcohorts, so the full
    # population support check is intentionally skipped for those calls.
    if (isTRUE(include_population)) {
      pop_b <- .recalc_population_weights_for_idx(
        month_df = month_df,
        pop_lookup = pop_lookup,
        idx = idx0,
        age_col = age_col,
        sex_col = sex_col,
        decile_col = decile_col,
        weight_cap = weight_cap
      )
      
      if (is.null(pop_b) || !nrow(pop_b)) {
        next
      }
      
      aligned <- ref_b %>%
        inner_join(pop_b, by = c("idx", "boot_pos")) %>%
        arrange(.data[["boot_pos"]]) %>%
        filter(
          is.finite(.data[["w_ref"]]), !is.na(.data[["w_ref"]]), .data[["w_ref"]] > 0,
          is.finite(.data[["w_pop"]]), !is.na(.data[["w_pop"]]), .data[["w_pop"]] > 0
        )
    } else {
      aligned <- ref_b %>%
        arrange(.data[["boot_pos"]]) %>%
        mutate(w_pop = NA_real_) %>%
        filter(
          is.finite(.data[["w_ref"]]), !is.na(.data[["w_ref"]]), .data[["w_ref"]] > 0
        )
    }
    
    if (!nrow(aligned)) {
      next
    }
    
    used <- used + 1L
    
    reps[[used]] <- list(
      idx = as.integer(aligned$idx),
      w_ref = as.numeric(aligned$w_ref),
      w_pop = as.numeric(aligned$w_pop),
      bootstrap_rep = b,
      n_sampled = length(idx0),
      n_used = nrow(aligned),
      has_population_weights = isTRUE(include_population)
    )
  }
  
  if (used < 1L) {
    return(vector("list", 0L))
  }
  
  reps[seq_len(used)]
}


# Logging wrapper used by run_models().
# This keeps the existing .safe_make_shared_reps() function name so the main
# runner does not need to change, but it no longer catches errors. If bootstrap
# construction fails unexpectedly, the pipeline now stops and shows the real
# error.
.safe_make_shared_reps <- function(month_df,
                                   ref_df,
                                   pop_lookup,
                                   boot_n,
                                   boot_seed,
                                   train_key,
                                   label,
                                   include_population = TRUE) {
  reps <- make_shared_bootstrap_reps(
    month_df = month_df,
    ref_df = ref_df,
    pop_lookup = pop_lookup,
    boot_n = boot_n,
    boot_seed = boot_seed,
    ref_month = train_key,
    age_col = "age_weight_raw",
    sex_col = "sexM",
    decile_col = "decile_weight_raw",
    ref_prefix = paste0("boot_w_ref_", gsub("[^A-Za-z0-9]+", "_", label)),
    weight_cap = 3.5,
    include_population = include_population
  )
  
  message(sprintf(
    "[run_models] Shared bootstrap %s: requested %d, successful %d%s",
    label,
    as.integer(boot_n),
    length(reps),
    ifelse(isTRUE(include_population), "", " (reference weights only)")
  ))
  
  reps
}
