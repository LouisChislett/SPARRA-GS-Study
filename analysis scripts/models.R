# =============================================================================
# models.R
#
# Author       : Louis Chislett
#
# Run the SPARRA v3 fixed odds-ratio model evaluation pipeline.
#
# This version always runs the full pipeline:
#   1. Read and standardise the input data.
#   2. Build the FE, LTC, YED, and Overall evaluation cohorts.
#   3. Score the fixed odds-ratio models for each evaluation month.
#   4. Create shared monthly bootstrap replicates.
#   5. Write metric outputs.
#   6. Write SHAP / exact contribution outputs.
#
# Compared with the original models.R, this version removes:
#   - outputs = c("metrics", "shap")
#   - selected prediction extracts
#   - long console summary prints
#   - tryCatch wrappers around reading and SHAP
#   - repeated FE/LTC/YED component-table code
#
# Assumptions:
#   - The fixed OR model objects and scoring helpers have already been sourced.
#   - make_sparra_v3_datasets() returns FE, LTC, YED, and Overall datasets.
#   - add_all_weights() has already created the required weight columns inside
#     make_sparra_v3_datasets().
#   - If something fails, the script should stop rather than silently continue.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(lubridate)
  library(fst)
})


# -----------------------------------------------------------------------------
# General utilities
# -----------------------------------------------------------------------------

# Write a data frame as UTF-8 CSV using consistent settings everywhere.
write_csv_utf8 <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE, fileEncoding = "UTF-8")
}


# Ensure every dataset has the month fields used by the monthly evaluation loop.
# make_sparra_v3_datasets() keeps `time`, but does not keep `year_month` or `date`,
# so this function reconstructs them when needed.
add_month_fields <- function(df) {
  if (!nrow(df)) return(df)
  
  df %>%
    mutate(
      year_month = if (!"year_month" %in% names(.)) {
        format(as_datetime(time, tz = "UTC"), "%Y-%m")
      } else year_month,
      date = if (!"date" %in% names(.)) {
        as.Date(paste0(year_month, "-01"))
      } else as.Date(date)
    ) %>%
    arrange(date) %>%
    mutate(across(where(is.character), ascii_clean))
}


# Fail early if a cohort is missing a column that the runner relies on later.
# This is intentionally stricter than silently creating fallback weights:
# if required model weights are missing, the input/cohort builder should be fixed.
check_required_cols <- function(df, cols, label) {
  missing_cols <- setdiff(cols, names(df))
  
  if (length(missing_cols)) {
    stop(sprintf(
      "%s is missing required column(s): %s",
      label,
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  invisible(df)
}


# Read the fst input once and do only the cleaning that is shared by all cohorts..
read_model_input <- function(in_file) {
  read_fst(in_file) %>%
    mutate(
      year_month = if (!"year_month" %in% names(.)) {
        format(as_datetime(time, tz = "UTC"), "%Y-%m")
      } else year_month,
      date = if (!"date" %in% names(.)) {
        as.Date(paste0(year_month, "-01"))
      } else as.Date(date),
      target = to01(target)
    ) %>%
    arrange(date) %>%
    mutate(across(where(is.character), ascii_clean))
}


# Build the four runner-ready cohorts.
# make_sparra_v3_datasets() does the substantive cohort work. This wrapper only:
#   - removes any remaining rows with missing target,
#   - restores year_month/date if needed,
#   - checks that expected weight columns exist.
prepare_cohorts <- function(df_all, pop_csv) {
  cohorts <- make_sparra_v3_datasets(df_all, pop_csv = pop_csv)
  
  cohorts <- lapply(cohorts, function(x) {
    x %>%
      filter(!is.na(target)) %>%
      add_month_fields()
  })
  
  check_required_cols(cohorts$FE,      c("year_month", "date", "target", "w_age_FE_trim"),  "FE")
  check_required_cols(cohorts$LTC,     c("year_month", "date", "target", "w_age_LTC_trim"), "LTC")
  check_required_cols(cohorts$YED,     c("year_month", "date", "target", "w_age_YED_trim"), "YED")
  check_required_cols(cohorts$Overall, c("year_month", "date", "target", "w_age_trim", "w_pop_trim"), "Overall")
  
  cohorts
}


# Store the small amount of metadata needed to run the same logic for FE/LTC/YED.
# The model name is also the scorer name, so score_model(mdl, dm) is used later.
make_model_specs <- function(cohorts) {
  list(
    FE  = list(df = cohorts$FE,  w_col = "w_age_FE_trim"),
    LTC = list(df = cohorts$LTC, w_col = "w_age_LTC_trim"),
    YED = list(df = cohorts$YED, w_col = "w_age_YED_trim")
  )
}


# Empty component table used when a model-specific cohort has no rows in a month.
# build_combined_on_overall() expects FE/LTC/YED tables with the same columns,
# even when one component is unavailable.
empty_component_tbl <- function() {
  tibble(id = integer(), y = integer(), p = double(), w = double())
}


# Convert a scored FE/LTC/YED monthly part into the standard component table used
# by build_combined_on_overall().
component_tbl <- function(part, w_col) {
  if (is.null(part)) return(empty_component_tbl())
  
  tibble(
    id = part$dm$id,
    y  = part$dm$target,
    p  = part$pred,
    w  = part$dm[[w_col]]
  )
}


# Calculate one metric row/table from a shared bootstrap replicate object.
# This function keeps the repeated metric_cis_from_shared_bootstrap() arguments
# in one place for FE/LTC/YED, COMBINED, and COMBINED_POP.
metric_row <- function(y, p, w_point, reps, weight, method, model, ym, weighting) {
  metric_cis_from_shared_bootstrap(
    y = y,
    p = p,
    w_point = w_point,
    reps = reps,
    weight = weight,
    thr = 0.5,
    conf.level = 0.95,
    bootstrap_method = method,
    return_boot = TRUE
  )$summary %>%
    mutate(model = model, year_month = ym, weighting = weighting)
}


# Produce row-level exact contributions and the corresponding weighted monthly
# SHAP summary for one model/month.
# The same monthly bootstrap replicates used for metrics are passed into the SHAP
# summary to keep uncertainty estimates aligned.
shap_month_outputs <- function(dm, mdl, w_col, ref_df, train_key, shap_ci_conf, reps, ym) {
  rows <- compute_exact_contributions(
    df_month = dm,
    mdl = mdl,
    w_col = w_col,
    ref_df = ref_df,
    ref_w_col = w_col
  ) %>%
    mutate(model = mdl)
  
  summary <- summarise_shap_weighted(
    rows,
    w_col = w_col,
    mdl = mdl,
    ref_df = ref_df,
    ref_month = train_key,
    conf = shap_ci_conf,
    shared_reps = setNames(list(reps), ym),
    weight = "w_ref"
  )
  
  list(rows = rows, summary = summary)
}


# Build the combined model output for one evaluation month.
# The combined model is evaluated on the Overall cohort. The helper maps the
# FE/LTC/YED predictions onto Overall using build_combined_on_overall().
make_combined_month <- function(overall_df, scored_parts, ym) {
  ov <- overall_df %>% filter(year_month == ym)
  
  build_combined_on_overall(
    overall_tbl = ov,
    fe_tbl  = component_tbl(scored_parts$FE,  "w_age_FE_trim"),
    ltc_tbl = component_tbl(scored_parts$LTC, "w_age_LTC_trim"),
    yed_tbl = component_tbl(scored_parts$YED, "w_age_YED_trim")
  )
}


# Keep only the fields needed to create combined-model bootstrap replicates.
# The combined replicates need both reference and population weights, so they are
# created from the Overall population rather than from FE/LTC/YED subcohorts.
make_combined_boot_df <- function(comb) {
  comb %>%
    transmute(
      id = id,
      time = time,
      target = y_overall,
      age_weight_raw = age_weight_raw,
      decile_weight_raw = decile_weight_raw,
      sexM = sexM,
      ..pred.. = p_combined
    )
}


# Finalise and write all metric outputs.
# Reference-weighted metrics include FE, LTC, YED, and COMBINED.
# Population-weighted metrics include only COMBINED_POP.
write_metric_outputs <- function(perf_rows_ref, perf_rows_pop, all_keys, out_dir) {
  out_dir_ref <- file.path(out_dir, "reference weighted")
  out_dir_pop <- file.path(out_dir, "population weighted")
  
  perf_df_ref <- bind_rows(perf_rows_ref) %>%
    left_join(all_keys, by = "year_month") %>%
    mutate(model = ascii_vec(model), weighting = ascii_vec(weighting)) %>%
    arrange(date, model)
  
  perf_df_pop <- bind_rows(perf_rows_pop) %>%
    left_join(all_keys, by = "year_month") %>%
    mutate(model = ascii_vec(model), weighting = ascii_vec(weighting)) %>%
    arrange(date, model)
  
  perf_df <- bind_rows(perf_df_ref, perf_df_pop) %>%
    arrange(date, weighting, model)
  
  write_csv_utf8(perf_df_ref, file.path(out_dir_ref, "performance_reference_weighted.csv"))
  write_csv_utf8(perf_df_pop, file.path(out_dir_pop, "performance_population_weighted.csv"))
  write_csv_utf8(perf_df,     file.path(out_dir,     "performance_all.csv"))
  
  list(
    performance = perf_df_ref,
    performance_population = perf_df_pop,
    performance_all = perf_df
  )
}


# Finalise and write all SHAP / exact contribution outputs.
# Each model gets:
#   - row-level contribution output,
#   - monthly weighted summary output,
#   - pre/post difference summary,
#   - pre/post bootstrap difference output.
write_shap_outputs <- function(shap_rows_by_model,
                               shap_summary_by_model,
                               shared_reps_by_model,
                               cohort_specs,
                               ref_dfs,
                               train_key,
                               pre_covid_key,
                               post_covid_key,
                               shap_ci_conf,
                               shap_dir) {
  shap_rows_all <- list()
  shap_summary_all <- list()
  diff_summary_by_model <- list()
  diff_boot_by_model <- list()
  
  for (mdl in names(cohort_specs)) {
    spec <- cohort_specs[[mdl]]
    ref_df <- ref_dfs[[mdl]]
    
    # If the reference month is unavailable for a model, fall back to all rows in
    # that model cohort. This preserves the original fallback behaviour.
    if (is.null(ref_df) || !nrow(ref_df)) {
      ref_df <- spec$df
    }
    
    shap_rows_mdl <- bind_rows(shap_rows_by_model[[mdl]])
    shap_summary_mdl <- bind_rows(shap_summary_by_model[[mdl]])
    
    diff_obj <- compute_shap_difference_weighted(
      shap_df = shap_rows_mdl,
      w_col = spec$w_col,
      mdl = mdl,
      ref_df = ref_df,
      ref_month = train_key,
      pre_month = pre_covid_key,
      post_month = post_covid_key,
      conf = shap_ci_conf,
      shared_reps = shared_reps_by_model[[mdl]],
      weight = "w_ref"
    )
    
    mdl_dir <- file.path(shap_dir, mdl)
    dir.create(mdl_dir, showWarnings = FALSE, recursive = TRUE)
    
    write_csv_utf8(shap_rows_mdl,    file.path(mdl_dir, sprintf("shap_rows_%s.csv", mdl)))
    write_csv_utf8(shap_summary_mdl, file.path(mdl_dir, sprintf("shap_summary_%s.csv", mdl)))
    
    write_csv_utf8(
      diff_obj$summary,
      file.path(
        mdl_dir,
        sprintf(
          "shap_diff_summary_%s_%s_minus_%s.csv",
          mdl,
          gsub("-", "_", post_covid_key),
          gsub("-", "_", pre_covid_key)
        )
      )
    )
    
    write_csv_utf8(
      diff_obj$boot,
      file.path(
        mdl_dir,
        sprintf(
          "shap_diff_bootstrap_%s_%s_minus_%s.csv",
          mdl,
          gsub("-", "_", post_covid_key),
          gsub("-", "_", pre_covid_key)
        )
      )
    )
    
    shap_rows_all[[mdl]] <- shap_rows_mdl
    shap_summary_all[[mdl]] <- shap_summary_mdl
    diff_summary_by_model[[mdl]] <- diff_obj$summary
    diff_boot_by_model[[mdl]] <- diff_obj$boot
    
    message(sprintf("[run_models] Saved SHAP tables for %s", mdl))
  }
  
  write_csv_utf8(
    bind_rows(shap_summary_all),
    file.path(shap_dir, "shap_summary_all_models.csv")
  )
  
  write_csv_utf8(
    bind_rows(diff_summary_by_model),
    file.path(
      shap_dir,
      sprintf(
        "shap_diff_summary_all_models_%s_minus_%s.csv",
        gsub("-", "_", post_covid_key),
        gsub("-", "_", pre_covid_key)
      )
    )
  )
  
  write_csv_utf8(
    bind_rows(diff_boot_by_model),
    file.path(
      shap_dir,
      sprintf(
        "shap_diff_bootstrap_all_models_%s_minus_%s.csv",
        gsub("-", "_", post_covid_key),
        gsub("-", "_", pre_covid_key)
      )
    )
  )
  
  list(
    rows = shap_rows_all,
    summary = shap_summary_all,
    diff_summary = diff_summary_by_model,
    diff_boot = diff_boot_by_model
  )
}


# -----------------------------------------------------------------------------
# Main runner
# -----------------------------------------------------------------------------

run_models <- function(in_file,
                       out_dir = "results",
                       train_key = "2012-01",
                       boot_n = 100,
                       boot_seed = 123,
                       pop_csv = POP_STANDARDISATION_CSV,
                       pre_covid_key = "2020-03",
                       post_covid_key = "2021-03",
                       shap_ci_conf = 0.95) {

  boot_n <- as.integer(boot_n)[1]
  
  # Create the output folder structure. Metrics and SHAP are always written.
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_dir_ref <- file.path(out_dir, "reference weighted")
  out_dir_pop <- file.path(out_dir, "population weighted")
  shap_dir <- file.path(out_dir, "shap")
  
  dir.create(out_dir_ref, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir_pop, showWarnings = FALSE, recursive = TRUE)
  dir.create(shap_dir,    showWarnings = FALSE, recursive = TRUE)
  
  # Read data, build cohorts, and load the population lookup used when
  # recalculating population weights in bootstrap replicates.
  df_all <- read_model_input(in_file)
  cohorts <- prepare_cohorts(df_all, pop_csv)
  pop_lookup <- get_population_standardisation_lookup_cached(pop_csv = pop_csv)
  
  # Evaluation months are all months after the fixed training/reference month.
  all_keys <- df_all %>%
    distinct(year_month, date) %>%
    arrange(date)
  
  eval_keys_tbl <- all_keys %>%
    filter(date > as.Date(paste0(train_key, "-01")))
  
  # Metadata used to apply identical logic to FE, LTC, and YED.
  cohort_specs <- make_model_specs(cohorts)
  
  # Reference-month datasets used for bootstrap weight recalculation and SHAP
  # reference comparisons.
  ref_dfs <- list(
    FE  = cohorts$FE  %>% filter(year_month == train_key),
    LTC = cohorts$LTC %>% filter(year_month == train_key),
    YED = cohorts$YED %>% filter(year_month == train_key),
    COMBINED = cohorts$Overall %>%
      filter(year_month == train_key) %>%
      distinct(id, year_month, .keep_all = TRUE) %>%
      transmute(
        id = id,
        time = time,
        target = target,
        age_weight_raw = age_weight_raw,
        decile_weight_raw = decile_weight_raw,
        sexM = sexM
      )
  )
  
  # Short method labels written to the metric output tables.
  metric_method_ref <- "shared_month_bootstrap_recalc_reference_weights"
  metric_method_pop <- "shared_month_bootstrap_recalc_population_weights"
  
  # Accumulators populated month-by-month and written once at the end.
  perf_rows_ref <- list()
  perf_rows_pop <- list()
  
  shap_rows_by_model <- list(FE = list(), LTC = list(), YED = list())
  shap_summary_by_model <- list(FE = list(), LTC = list(), YED = list())
  shared_reps_by_model <- list(FE = list(), LTC = list(), YED = list())
  
  message(sprintf(
    "[run_models] Running metrics and SHAP | boot_n:%d boot_seed:%d",
    boot_n,
    as.integer(boot_seed)
  ))
  
  # ---------------------------------------------------------------------------
  # Monthly evaluation loop
  # ---------------------------------------------------------------------------
  for (j in seq_len(nrow(eval_keys_tbl))) {
    ym <- eval_keys_tbl$year_month[[j]]
    
    # Give each month a deterministic seed block so that each model/month gets a
    # reproducible but distinct bootstrap path.
    month_seed_base <- as.integer(boot_seed + j * 10000L)
    
    message(sprintf("[run_models] Month %s", ym))
    
    # Stores this month's component model predictions. These are needed to build
    # the combined model after FE/LTC/YED have been scored.
    scored_parts <- list()
    
    # Score FE, LTC, and YED for this month.
    for (i in seq_along(cohort_specs)) {
      mdl <- names(cohort_specs)[[i]]
      spec <- cohort_specs[[mdl]]
      dm <- spec$df %>% filter(year_month == ym)
      
      if (!nrow(dm)) next
      
      # 'Fit' the OR model.
      pred <- score_model(mdl, dm)
      
      # Build one shared bootstrap path for this model/month. The same replicate
      # object is used for metrics and SHAP summaries.
      reps <- .safe_make_shared_reps(
        month_df = dm %>% mutate(..pred.. = pred),
        ref_df = ref_dfs[[mdl]],
        pop_lookup = pop_lookup,
        boot_n = boot_n,
        boot_seed = month_seed_base + i,
        train_key = train_key,
        label = paste(mdl, ym, sep = "_"),
        include_population = FALSE
      )
      
      # this is used by SHAP later
      shared_reps_by_model[[mdl]][[ym]] <- reps
      scored_parts[[mdl]] <- list(dm = dm, pred = pred)
      
      # Reference-weighted metrics for the model-specific cohort.
      perf_rows_ref[[length(perf_rows_ref) + 1L]] <- metric_row(
        y = dm$target,
        p = pred,
        w_point = dm[[spec$w_col]],
        reps = reps,
        weight = "w_ref",
        method = metric_method_ref,
        model = mdl,
        ym = ym,
        weighting = "reference weighted"
      )
      
      # Row-level and summary SHAP outputs for the model-specific cohort.
      shap_obj <- shap_month_outputs(
        dm = dm,
        mdl = mdl,
        w_col = spec$w_col,
        ref_df = ref_dfs[[mdl]],
        train_key = train_key,
        shap_ci_conf = shap_ci_conf,
        reps = reps,
        ym = ym
      )
      
      shap_rows_by_model[[mdl]][[length(shap_rows_by_model[[mdl]]) + 1L]] <- shap_obj$rows
      shap_summary_by_model[[mdl]][[length(shap_summary_by_model[[mdl]]) + 1L]] <- shap_obj$summary
    }
    
    # Build the combined model on the Overall population after component scoring.
    comb <- make_combined_month(
      overall_df = cohorts$Overall,
      scored_parts = scored_parts,
      ym = ym
    )
    
    # Combined metrics need both reference and population bootstrap weights.
    reps_combined <- .safe_make_shared_reps(
      month_df = make_combined_boot_df(comb),
      ref_df = ref_dfs$COMBINED,
      pop_lookup = pop_lookup,
      boot_n = boot_n,
      boot_seed = month_seed_base + 100L,
      train_key = train_key,
      label = paste("COMBINED", ym, sep = "_"),
      include_population = TRUE
    )
    
    # Reference-weighted combined metrics.
    perf_rows_ref[[length(perf_rows_ref) + 1L]] <- metric_row(
      y = comb$y_overall,
      p = comb$p_combined,
      w_point = comb$w_overall,
      reps = reps_combined,
      weight = "w_ref",
      method = metric_method_ref,
      model = "COMBINED",
      ym = ym,
      weighting = "reference weighted"
    )
    
    # Population-weighted combined metrics.
    perf_rows_pop[[length(perf_rows_pop) + 1L]] <- metric_row(
      y = comb$y_overall,
      p = comb$p_combined,
      w_point = comb$w_pop_overall,
      reps = reps_combined,
      weight = "w_pop",
      method = metric_method_pop,
      model = "COMBINED_POP",
      ym = ym,
      weighting = "population weighted"
    )
  }
  
  # Write metric outputs after all months have been evaluated.
  metrics_out <- write_metric_outputs(
    perf_rows_ref = perf_rows_ref,
    perf_rows_pop = perf_rows_pop,
    all_keys = all_keys,
    out_dir = out_dir
  )
  
  # Write SHAP outputs after all monthly row-level and summary pieces exist.
  shap_out <- write_shap_outputs(
    shap_rows_by_model = shap_rows_by_model,
    shap_summary_by_model = shap_summary_by_model,
    shared_reps_by_model = shared_reps_by_model,
    cohort_specs = cohort_specs,
    ref_dfs = ref_dfs,
    train_key = train_key,
    pre_covid_key = pre_covid_key,
    post_covid_key = post_covid_key,
    shap_ci_conf = shap_ci_conf,
    shap_dir = shap_dir
  )
  
  # Return the final in-memory outputs for interactive use/tests. The main side
  # effect of the function is still the CSV files written to out_dir.
  invisible(list(
    metrics = metrics_out,
    shap = shap_out
  ))
}
