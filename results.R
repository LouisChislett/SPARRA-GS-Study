# =============================================================================
# results.R
#
# Author       : Louis Chislett
#
# Purpose :
#   Run the SPARRA analysis pipeline:
#
#     (1) EXPLORATORY PLOTS
#         - sources exploratoryPlots.R
#         - produces descriptive overview plots for the unweighted,
#           reference-month weighted, and population-weighted datasets
#
#     (2) MODEL EVALUATION / SHAP TABLES
#         - sources models.R
#         - loads ors.R when needed
#         - computes model performance metrics and exact contribution / SHAP
#           summary tables
#         - writes tabular outputs to the results folder
#
#     (3) MODEL PLOTS
#         - sources modelPlots.R
#         - rebuilds model and SHAP plots from the saved CSV outputs in the
#           results folder
#
#     (4) WEIGHTED DATA DRIFT
#         - sources weightedDataDrift.R
#         - compares two selected months and saves drift outputs
# =============================================================================

suppressPackageStartupMessages({
  library(fst)
  library(dplyr)
  library(lubridate)
})

# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
cfg <- list(
  in_file = "PATH/cleanData/covariate_matrix_monthly_with_target.fst",

  exploratory_out_dir    = "exploratory plots",
  model_out_dir          = "model plots",
  results_dir            = "results",
  drift_folder           = "data drift",
  temporal_refit_out_dir = "temporal refit comparison",

  # Folder containing reusable analysis scripts
  scripts_dir = "analysis scripts",

  seed = 123,

  drift_month1 = "2020-03-01",
  drift_month2 = "2021-03-01",

  train_key = "2012-01",
  topK      = 10,

  # Pre/post comparison months used by models.R and modelPlots.R
  pre_covid_key  = "2020-03",
  post_covid_key = "2021-03",

  temporal_refit_test_frac = 0.50,
  temporal_refit_include_train_year_eval = TRUE,

  run_exploratory    = TRUE,
  run_drift          = TRUE,
  run_model_eval     = TRUE,
  run_model_plots    = TRUE,
  run_shap           = TRUE,
  run_temporal_refit = TRUE,

  source_variable_checking = TRUE,

  # Model evaluation / exact contribution settings
  boot_n       = 100,
  boot_seed    = 123,
  shap_ci_conf = 0.95,
  results_csv  = "results.csv",

  cohort_builder_file = "cohortBuilder.R",
  exploratory_file    = "exploratoryPlots.R",
  weighted_drift_file = "weightedDataDrift.R",
  models_file         = "models.R",
  model_plots_file    = "modelPlots.R",
  ors_file            = "ors.R",
  metrics_file        = "metrics.R",
  shap_file           = "shap.R",
  utils_file          = "utils.R",
  bootstrap_file      = "bootstrapHelpers.R",
  models_scoring_file = "modelScoring.R",
  standardisation_file= "standardisationHelpers.R"
)

# -----------------------------------------------------------------------------
# Source project code
# -----------------------------------------------------------------------------
source(file.path(cfg$scripts_dir, cfg$cohort_builder_file))
source(file.path(cfg$scripts_dir, cfg$metrics_file))
source(file.path(cfg$scripts_dir, cfg$shap_file))
source(file.path(cfg$scripts_dir, cfg$utils_file))
source(file.path(cfg$scripts_dir, cfg$bootstrap_file))
source(file.path(cfg$scripts_dir, cfg$models_scoring_file))
source(file.path(cfg$scripts_dir, cfg$standardisation_file))
source(file.path(cfg$scripts_dir, cfg$exploratory_file))
source(file.path(cfg$scripts_dir, cfg$weighted_drift_file))
source(file.path(cfg$scripts_dir, cfg$models_file))
source(file.path(cfg$scripts_dir, cfg$model_plots_file))
source(file.path(cfg$scripts_dir, cfg$ors_file))

# -----------------------------------------------------------------------------
# STEP 1: Exploratory plots
# -----------------------------------------------------------------------------
message("[STEP 1] Exploratory Plots")
dir.create(cfg$exploratory_out_dir, recursive = TRUE, showWarnings = FALSE)

res_exploratory <- makeExploratoryPlots(
  in_file        = cfg$in_file,
  out_dir        = cfg$exploratory_out_dir,
  seed           = cfg$seed,
  weight_col     = "w_age_trim",
  weight_col_FE  = "w_age_FE_trim",
  weight_col_LTC = "w_age_LTC_trim",
  weight_col_YED = "w_age_YED_trim",
  weight_col_pop = "w_pop_trim"
)

message("[STEP 1] Done")


# -----------------------------------------------------------------------------
# STEP 2: Model evaluation and SHAP
# -----------------------------------------------------------------------------
message("[STEP 2] Model evaluation and SHAP tables")
dir.create(cfg$results_dir, recursive = TRUE, showWarnings = FALSE)

res_models <- run_models(
  in_file        = cfg$in_file,
  out_dir        = cfg$results_dir,
  train_key      = cfg$train_key,
  boot_n         = cfg$boot_n,
  boot_seed      = cfg$boot_seed,
  pop_csv        = POP_STANDARDISATION_CSV,
  pre_covid_key  = cfg$pre_covid_key,
  post_covid_key = cfg$post_covid_key,
  shap_ci_conf   = cfg$shap_ci_conf
)
  
message("[STEP 2] Done")


# -----------------------------------------------------------------------------
# STEP 2B: Model plots from saved CSVs
# -----------------------------------------------------------------------------
message("[STEP 2C] Model plots from saved CSVs")
dir.create(cfg$model_out_dir, recursive = TRUE, showWarnings = FALSE)

res_model_plots <- run_model_plots(
  results_dir    = cfg$results_dir,
  plot_dir       = cfg$model_out_dir,
  train_key      = cfg$train_key,
  shap_ci_conf   = cfg$shap_ci_conf,
  pre_covid_key  = cfg$pre_covid_key,
  post_covid_key = cfg$post_covid_key
)

message("[STEP 2B] Done")

# -----------------------------------------------------------------------------
# STEP 3: Weighted data drift analysis
# -----------------------------------------------------------------------------
message("[STEP 3] Weighted data drift analysis")
dir.create(cfg$drift_folder, recursive = TRUE, showWarnings = FALSE)

res_drift <- weightedDataDrift(
  in_file = cfg$in_file,
  month1 = cfg$drift_month1,
  month2 = cfg$drift_month2,
  out_dir = cfg$drift_folder,
  include_overall = FALSE,
  make_plots = TRUE,
  summary_csv_name = "resultsdataDrift.csv"
)

message("[STEP 3] Done")

