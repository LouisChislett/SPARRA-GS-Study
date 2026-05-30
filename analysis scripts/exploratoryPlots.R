# =============================================================================
# exploratoryPlots.R
#
# Author       : Louis Chislett
#
# Purpose     :
#   Generate exploratory plots for SPARRA v3 GS study
#
#   CURRENT PLOTTING CONVENTION:
#     - cohorts_w   = RAW data with analysis weights retained
#     - cohorts_ns  = RAW data used for UNWEIGHTED plots, with analysis-weight
#                     columns removed before plotting
#
#   Weighting columns expected on weighted datasets:
#     - Overall age-weighted      : w_age_trim
#     - FE age-weighted           : w_age_FE_trim
#     - LTC age-weighted          : w_age_LTC_trim
#     - YED age-weighted          : w_age_YED_trim
#     - Overall population-weight : w_pop_trim
#
# =============================================================================

makeExploratoryPlots <- function(
    in_file = "PATH/cleanData/covariate_matrix_monthly_with_target.fst",
    out_dir = "exploratory plots",
    pop_csv = "scotland_midyear_pop_age_sex_simd_decile_2012_2022_no_persons.csv",
    seed   = 123,
    weight_col = "w_age_trim",
    weight_col_FE  = "w_age_FE_trim",
    weight_col_LTC = "w_age_LTC_trim",
    weight_col_YED = "w_age_YED_trim",
    weight_col_pop = "w_pop_trim"
) {
  suppressPackageStartupMessages({
    library(fst)
    library(dplyr)
    library(lubridate)
    library(ggplot2)
    library(stringr)
    library(scales)
    library(readr)
    library(tidyr)
  })
  set.seed(seed)
  
  # ---- source dependencies (project functions) ----
  if (file.exists("variableChecking.R")) source("variableChecking.R")
  
  if (!exists("make_sparra_v3_datasets_raw")) {
    if (file.exists("make_sparra_v3_datasets_raw.R")) source("make_sparra_v3_datasets_raw.R")
    else stop("make_sparra_v3_datasets_raw() not found.")
  }
  # ---- output dirs ----
  dirs <- list(
    overview_w     = file.path(out_dir, "overview", "weighted"),
    overview_ns    = file.path(out_dir, "overview", "non_stratified"),
    overview_pop   = file.path(out_dir, "overview", "population_weighted"),
    variables      = file.path(out_dir, "variables"),
    weight_diag    = file.path(out_dir, "weights diagnostics"),
    age_w          = file.path(out_dir, "age", "weighted"),
    age_ns         = file.path(out_dir, "age", "non_stratified"),
    age_pop        = file.path(out_dir, "age", "population_weighted"),
    tr_w_g         = file.path(out_dir, "target_reasons", "weighted", "grouped"),
    tr_ns_g        = file.path(out_dir, "target_reasons", "non_stratified", "grouped"),
    tr_pop_g       = file.path(out_dir, "target_reasons", "population_weighted", "grouped"),
    decile         = file.path(out_dir, "decile")
  )
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)
  
  # ---- helpers ----
  to01 <- function(x) {
    if (is.logical(x)) return(as.integer(x))
    if (is.numeric(x)) {
      ux <- unique(na.omit(x))
      if (all(ux %in% c(0, 1))) return(as.integer(x))
      stop("Numeric target has non-binary values; refusing to coerce.")
    }
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) {
      tx <- toupper(trimws(x))
      return(ifelse(tx %in% c("1","TRUE","T","YES","Y"), 1L,
                    ifelse(tx %in% c("0","FALSE","F","NO","N"), 0L, NA_integer_)))
    }
    stop("Unsupported target type: ", paste(class(x), collapse = ", "))
  }
  
  add_dates <- function(df) {
    df %>%
      mutate(
        year_month = if (!"year_month" %in% names(.)) format(as_datetime(time, tz = "UTC"), "%Y-%m") else year_month,
        date       = if (!"date" %in% names(.)) as.Date(paste0(year_month, "-01")) else as.Date(date)
      ) %>%
      arrange(date)
  }
  
  safe_filename <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)
  
  parse_population_age_to_numeric <- function(x) {
    y <- trimws(as.character(x))
    y_low <- tolower(y)
    dplyr::case_when(
      y_low %in% c("total") ~ NA_real_,
      grepl("^\\d+\\+$", y) ~ as.numeric(sub("\\+$", "", y)),
      grepl("^\\d+$", y)    ~ as.numeric(y),
      TRUE                  ~ suppressWarnings(as.numeric(y))
    )
  }
  
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
  
  normalise_within_source <- function(df, value_col = "value") {
    df %>%
      group_by(source) %>%
      mutate(
        total_value = sum(.data[[value_col]], na.rm = TRUE),
        !!value_col := ifelse(total_value > 0, .data[[value_col]] / total_value, NA_real_)
      ) %>%
      ungroup() %>%
      select(-total_value)
  }
  
  make_compare_labels <- function(ref_month) {
    ref_year <- year(ref_month)
    ref_label <- format(ref_month, "%Y-%m")
    c(
      unweighted = paste0(ref_label, " Unweighted"),
      weighted   = paste0(ref_label, " Weighted"),
      population = paste0("Mid-Year Population ", ref_year)
    )
  }
  
  make_compare_palette <- function(ref_month) {
    labels <- unname(make_compare_labels(ref_month))
    vals <- c("#0072B2", "#D55E00", "#009E73")
    names(vals) <- labels
    vals
  }
  
  prop_axis_upper <- function(x, expand_mult = 1.06, min_pad = 0.01, max_upper = 1) {
    x <- x[is.finite(x) & !is.na(x)]
    if (length(x) == 0) return(max_upper)
    x_max <- max(x)
    if (x_max <= 0) return(0.05)
    upper <- x_max * expand_mult
    upper <- max(upper, x_max + min_pad)
    upper <- min(max_upper, upper)
    upper
  }
  
  # ---- weight resolution ----
  weight_map <- c(
    Overall = weight_col,
    FE      = weight_col_FE,
    LTC     = weight_col_LTC,
    YED     = weight_col_YED
  )
  
  resolve_weight_col <- function(cohort = NULL) {
    if (!is.null(cohort) && cohort %in% names(weight_map)) return(unname(weight_map[[cohort]]))
    weight_col
  }
  
  # All final analysis-weight columns expected on the weighted/raw datasets.
  # These are handled separately in the weights diagnostics below and are
  # excluded from the ordinary weighted-mean variable plots.
  weight_cols_all <- unique(c(unname(weight_map), weight_col_pop))
  
  is_analysis_weight_var <- function(v) {
    v %in% weight_cols_all ||
      grepl("^w_(age|age_FE|age_LTC|age_YED|pop)_(raw|stab|trim)$", v)
  }
  
  strip_analysis_weights <- function(df) {
    # Raw datasets include the same analysis-weight columns used for weighted plots.
    # For all formerly unstratified/unweighted plots, remove those columns so helper
    # functions fall back to equal weights via get_w_from_col(...).
    dplyr::select(df, -dplyr::any_of(weight_cols_all))
  }
  
  get_w_from_col <- function(df, wc) {
    if (!is.null(wc) && wc %in% names(df)) {
      w <- suppressWarnings(as.numeric(df[[wc]]))
      w[!is.finite(w)] <- NA_real_
      w[is.na(w)] <- 0
      w[w < 0] <- 0
      return(w)
    }
    rep(1, nrow(df))
  }
  
  get_w <- function(df, cohort = NULL) {
    get_w_from_col(df, resolve_weight_col(cohort))
  }
  
  # ---- plotting defaults ----
  theme_base <- theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 20, face = "bold", colour = "black"),
      axis.title.y = element_text(size = 20, face = "bold", colour = "black"),
      axis.text.x  = element_text(size = 16, colour = "black", angle = 90, vjust = 1, hjust = 1),
      axis.text.y  = element_text(size = 16, colour = "black"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.title    = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  
  scale_date_6m <- scale_x_date(
    date_breaks = "6 months",
    date_labels = "%Y-%b",
    expand = expansion(mult = c(0.005, 0.005))
  )
  
  save_plot <- function(path, p, w = 10, h = 5.5) {
    ggsave(filename = path, plot = p, width = w, height = h, dpi = 300, bg = "white")
  }
  
  # ---- weighted summaries ----
  weighted_prev <- function(target_vec, w) {
    ok <- !is.na(target_vec)
    denom <- sum(w[ok], na.rm = TRUE)
    if (denom <= 0) return(NA_real_)
    num <- sum(w[ok & target_vec == 1L], na.rm = TRUE)
    num / denom
  }
  
  ess_from_w <- function(w) {
    w <- w[is.finite(w) & !is.na(w)]
    if (length(w) == 0) return(NA_real_)
    (sum(w)^2) / sum(w^2)
  }
  
  monthly_prev <- function(df, cohort) {
    w <- get_w(df, cohort = cohort)
    df %>%
      mutate(w = w) %>%
      group_by(date) %>%
      summarise(
        n_total = n(),
        n_label = sum(!is.na(target)),
        prop_pos = weighted_prev(target, w),
        ess = ess_from_w(w[!is.na(target)]),
        .groups = "drop"
      ) %>%
      mutate(cohort = cohort)
  }
  
  monthly_prev_by_weightcol <- function(df, cohort_label, wc = NULL) {
    w <- get_w_from_col(df, wc)
    df %>%
      mutate(w = w) %>%
      group_by(date) %>%
      summarise(
        n_total = n(),
        n_label = sum(!is.na(target)),
        prop_pos = weighted_prev(target, w),
        ess = ess_from_w(w[!is.na(target)]),
        .groups = "drop"
      ) %>%
      mutate(cohort = cohort_label)
  }
  
  monthly_size <- function(df, cohort) {
    df %>%
      group_by(date) %>%
      summarise(n = n_distinct(id), .groups = "drop") %>%
      mutate(cohort = cohort)
  }
  
  is_binary_vec <- function(x) {
    ux <- unique(x[!is.na(x)])
    length(ux) > 0 && length(ux) <= 2 && all(ux %in% c(0, 1))
  }
  
  summarise_numeric <- function(df, v, cohort) {
    w <- get_w(df, cohort = cohort)
    df %>%
      mutate(w = w) %>%
      group_by(date) %>%
      summarise(
        mean_val = weighted.mean(.data[[v]], w, na.rm = TRUE),
        ess = ess_from_w(w[!is.na(.data[[v]])]),
        .groups = "drop"
      ) %>%
      mutate(cohort = cohort, var = v)
  }
  
  # ---- palettes ----
  # Colourblind-safe palette with explicit unnamed hex values to avoid ggplot
  # scale matching issues.
  okabe_ito <- c(
    black          = "#000000",
    orange         = "#E69F00",
    sky_blue       = "#56B4E9",
    bluish_green   = "#009E73",
    yellow         = "#F0E442",
    blue           = "#0072B2",
    vermillion     = "#D55E00",
    reddish_purple = "#CC79A7"
  )
  
  cohort_palette <- c(
    Overall = "#0072B2",
    FE      = "#D55E00",
    LTC     = "#009E73",
    YED     = "#CC79A7"
  )
  
  sub_palette <- c(
    FE  = "#D55E00",
    LTC = "#009E73",
    YED = "#CC79A7"
  )
  
  summarise_mean_weight <- function(df, cohort, wc) {
    if (is.null(df) || !nrow(df) || !wc %in% names(df)) {
      return(tibble(
        date = as.Date(character()),
        cohort = character(),
        weight_col = character(),
        n = integer(),
        n_non_missing = integer(),
        mean_weight = numeric()
      ))
    }
    
    df %>%
      transmute(
        date = as.Date(date),
        `..weight_value..` = suppressWarnings(as.numeric(.data[[wc]]))
      ) %>%
      group_by(date) %>%
      summarise(
        n = n(),
        n_non_missing = sum(is.finite(`..weight_value..`) & !is.na(`..weight_value..`)),
        mean_weight = {
          ww <- `..weight_value..`[is.finite(`..weight_value..`) & !is.na(`..weight_value..`)]
          if (length(ww) > 0) mean(ww) else NA_real_
        },
        .groups = "drop"
      ) %>%
      filter(n_non_missing > 0) %>%
      mutate(cohort = cohort, weight_col = wc) %>%
      select(date, cohort, weight_col, n, n_non_missing, mean_weight)
  }
  
  save_weight_diagnostic_plots <- function(cohort_list, out_dir_weight_diag) {
    dir.create(out_dir_weight_diag, recursive = TRUE, showWarnings = FALSE)
    
    # Start clean every run so old multi-cohort diagnostics cannot be mistaken
    # for the newly produced single-cohort checks.
    unlink(Sys.glob(file.path(out_dir_weight_diag, "*.png")))
    unlink(file.path(out_dir_weight_diag, "mean_weights_by_month_cohort.csv"))
    
    # IMPORTANT:
    # Each weight is checked ONLY in the cohort it was normalised for.
    # This deliberately avoids four-line plots such as w_age_FE_trim shown for
    # Overall/FE/LTC/YED. The purpose here is just: does the relevant mean stay
    # close to 1 over time?
    weight_diag_targets <- tibble::tibble(
      weight_col = c(weight_col, weight_col_FE, weight_col_LTC, weight_col_YED, weight_col_pop),
      cohort     = c("Overall",  "FE",          "LTC",          "YED",          "Overall"),
      file_stub  = c("w_age_trim_Overall",
                     "w_age_FE_trim_FE",
                     "w_age_LTC_trim_LTC",
                     "w_age_YED_trim_YED",
                     "w_pop_trim_Overall")
    ) %>%
      filter(!is.na(weight_col), weight_col != "") %>%
      distinct(weight_col, cohort, file_stub)
    
    diag_parts <- list()
    
    for (i in seq_len(nrow(weight_diag_targets))) {
      wc <- weight_diag_targets$weight_col[[i]]
      nm <- weight_diag_targets$cohort[[i]]
      stub <- weight_diag_targets$file_stub[[i]]
      
      if (is.null(cohort_list[[nm]]) || !wc %in% names(cohort_list[[nm]])) next
      
      dd <- summarise_mean_weight(cohort_list[[nm]], cohort = nm, wc = wc)
      if (!nrow(dd)) next
      
      diag_parts[[length(diag_parts) + 1L]] <- dd
      line_col <- if (nm %in% names(cohort_palette)) unname(cohort_palette[[nm]]) else "black"
      
      p <- ggplot(dd, aes(x = date, y = mean_weight)) +
        geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
        geom_line(colour = line_col, linewidth = 0.7, show.legend = FALSE) +
        geom_point(colour = line_col, size = 1.2, show.legend = FALSE) +
        scale_date_6m +
        scale_y_continuous(labels = label_number(accuracy = 0.001)) +
        labs(
          x = "Month",
          y = paste0("Mean ", wc),
          caption = paste0("Cohort checked: ", nm, "; dashed line = 1")
        ) +
        theme_base +
        theme(
          legend.position = "none",
          plot.caption = element_text(size = 10, colour = "black")
        )
      
      save_plot(
        file.path(out_dir_weight_diag, paste0("mean_", safe_filename(stub), ".png")),
        p, w = 10, h = 5.5
      )
    }
    
    diag <- bind_rows(diag_parts) %>% arrange(weight_col, cohort, date)
    if (nrow(diag)) {
      readr::write_csv(diag, file.path(out_dir_weight_diag, "mean_weights_by_month_cohort.csv"))
    } else {
      warning("No target-cohort analysis weights found; skipping weights diagnostics.", call. = FALSE)
    }
    
    invisible(diag)
  }
  
  grouped_levels <- c(
    "Certain infectious and parasitic diseases",
    "Neoplasms",
    "Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism",
    "Endocrine, nutritional and metabolic diseases",
    "Mental and behavioural disorders",
    "Diseases of the nervous system",
    "Diseases of the eye and adnexa",
    "Diseases of the ear and mastoid process",
    "Diseases of the circulatory system",
    "Diseases of the respiratory system",
    "Diseases of the digestive system",
    "Diseases of the skin and subcutaneous tissue",
    "Diseases of the musculoskeletal system and connective tissue",
    "Diseases of the genitourinary system",
    "Pregnancy, childbirth and the puerperium",
    "Certain conditions originating in the perinatal period",
    "Congenital malformations, deformations and chromosomal abnormalities",
    "Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified",
    "Injury, poisoning and certain other consequences of external causes",
    "External causes of morbidity and mortality",
    "Factors influencing health status and contact with health services",
    "Codes for special purposes",
    "OTHER_OR_UNMAPPED"
  )
  
  grouped_labels <- c(
    "Certain infectious and parasitic diseases" =
      "Infectious/parasitic",
    "Neoplasms" =
      "Neoplasms",
    "Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism" =
      "Blood/immune",
    "Endocrine, nutritional and metabolic diseases" =
      "Endocrine/metabolic",
    "Mental and behavioural disorders" =
      "Mental/behavioural",
    "Diseases of the nervous system" =
      "Nervous system",
    "Diseases of the eye and adnexa" =
      "Eye/adnexa",
    "Diseases of the ear and mastoid process" =
      "Ear/mastoid",
    "Diseases of the circulatory system" =
      "Circulatory",
    "Diseases of the respiratory system" =
      "Respiratory",
    "Diseases of the digestive system" =
      "Digestive",
    "Diseases of the skin and subcutaneous tissue" =
      "Skin/subcutaneous",
    "Diseases of the musculoskeletal system and connective tissue" =
      "Musculoskeletal",
    "Diseases of the genitourinary system" =
      "Genitourinary",
    "Pregnancy, childbirth and the puerperium" =
      "Pregnancy/childbirth",
    "Certain conditions originating in the perinatal period" =
      "Perinatal",
    "Congenital malformations, deformations and chromosomal abnormalities" =
      "Congenital",
    "Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified" =
      "Symptoms/signs",
    "Injury, poisoning and certain other consequences of external causes" =
      "Injury/poisoning",
    "External causes of morbidity and mortality" =
      "External causes",
    "Factors influencing health status and contact with health services" =
      "Health status/contact",
    "Codes for special purposes" =
      "Special purposes",
    "OTHER_OR_UNMAPPED" =
      "Other/unmapped"
  )
  
  grouped_palette <- c(
    "Certain infectious and parasitic diseases" =
      "#1b9e77",
    "Neoplasms" =
      "#d95f02",
    "Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism" =
      "#7570b3",
    "Endocrine, nutritional and metabolic diseases" =
      "#e7298a",
    "Mental and behavioural disorders" =
      "#66a61e",
    "Diseases of the nervous system" =
      "#e6ab02",
    "Diseases of the eye and adnexa" =
      "#a6761d",
    "Diseases of the ear and mastoid process" =
      "#666666",
    "Diseases of the circulatory system" =
      "#1f78b4",
    "Diseases of the respiratory system" =
      "#33a02c",
    "Diseases of the digestive system" =
      "#fb9a99",
    "Diseases of the skin and subcutaneous tissue" =
      "#fdbf6f",
    "Diseases of the musculoskeletal system and connective tissue" =
      "#cab2d6",
    "Diseases of the genitourinary system" =
      "#6a3d9a",
    "Pregnancy, childbirth and the puerperium" =
      "#ff7f00",
    "Certain conditions originating in the perinatal period" =
      "#b15928",
    "Congenital malformations, deformations and chromosomal abnormalities" =
      "#8dd3c7",
    "Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified" =
      "#ffffb3",
    "Injury, poisoning and certain other consequences of external causes" =
      "#fb8072",
    "External causes of morbidity and mortality" =
      "#80b1d3",
    "Factors influencing health status and contact with health services" =
      "#b3de69",
    "Codes for special purposes" =
      "#fccde5",
    "OTHER_OR_UNMAPPED" =
      "#bdbdbd"
  )
  
  # ---- load data (restrict window) ----
  df_all <- read_fst(in_file) %>%
    mutate(
      year_month = if (!"year_month" %in% names(.)) format(as_datetime(time, tz = "UTC"), "%Y-%m") else year_month,
      date       = if (!"date" %in% names(.)) as.Date(paste0(year_month, "-01")) else as.Date(date),
      target     = to01(target)
    ) %>%
    arrange(date) %>%
    filter(date >= as.Date("2012-01-01"), date <= as.Date("2022-10-01"))
  
  # ---- load population data for comparison plots ----
  if (!file.exists(pop_csv)) {
    stop("Population CSV not found: ", pop_csv)
  }
  
  pop_raw <- readr::read_csv(pop_csv, show_col_types = FALSE)
  
  needed_pop_cols <- c("year", "sex", "simd_decile", "age", "population")
  missing_pop_cols <- setdiff(needed_pop_cols, names(pop_raw))
  if (length(missing_pop_cols) > 0) {
    stop("Population CSV missing required columns: ", paste(missing_pop_cols, collapse = ", "))
  }
  
  pop_clean <- pop_raw %>%
    mutate(
      year = suppressWarnings(as.integer(year)),
      simd_decile = suppressWarnings(as.integer(simd_decile)),
      age_num = parse_population_age_to_numeric(age),
      population = suppressWarnings(as.numeric(population)),
      sex_std = standardise_sex_from_population(sex)
    ) %>%
    filter(
      !is.na(year),
      !is.na(simd_decile),
      !is.na(age_num),
      !is.na(population)
    )
  
  # ---- build cohorts: weighted + raw-for-unweighted ----
  build_sets <- function(builder_fun) {
    coh <- builder_fun(df_all)
    FE      <- add_dates(coh$FE      %>% filter(!is.na(target)))
    LTC     <- add_dates(coh$LTC     %>% filter(!is.na(target)))
    YED     <- add_dates(coh$YED     %>% filter(!is.na(target)))
    Overall <- add_dates(coh$Overall %>% filter(!is.na(target)))
    
    list(FE = FE, LTC = LTC, YED = YED, Overall = Overall)
  }
  
  cohorts_w  <- build_sets(make_sparra_v3_datasets_raw)
  cohorts_ns <- lapply(build_sets(make_sparra_v3_datasets_raw), strip_analysis_weights)
  
  weight_diagnostics <- save_weight_diagnostic_plots(cohorts_w, dirs$weight_diag)
  
  # =============================================================================
  # 1) SIMD decile comparison:
  #    2012-01 unweighted vs 2012-01 weighted vs mid-year population 2012
  #    SHOWS PROPORTIONS
  # =============================================================================
  plot_simd_population_compare <- function(df_w_overall, df_ns_overall, pop_df, out_path) {
    if (!("decile" %in% names(df_w_overall)) || !("decile" %in% names(df_ns_overall))) {
      warning("decile not found; skipping SIMD population comparison: ", out_path)
      return(invisible(NULL))
    }
    
    ref_month <- as.Date("2012-01-01")
    ref_year  <- 2012L
    compare_labels <- make_compare_labels(ref_month)
    compare_palette <- make_compare_palette(ref_month)
    
    dat_uw <- df_ns_overall %>%
      filter(date == ref_month, !is.na(decile)) %>%
      mutate(
        decile_chr = as.character(decile),
        source = compare_labels["unweighted"],
        value = 1
      ) %>%
      group_by(source, decile_chr) %>%
      summarise(value = sum(value), .groups = "drop")
    
    dat_w <- df_w_overall %>%
      filter(date == ref_month, !is.na(decile)) %>%
      mutate(
        decile_chr = as.character(decile),
        source = compare_labels["weighted"],
        w = get_w_from_col(., weight_col_pop)
      ) %>%
      group_by(source, decile_chr) %>%
      summarise(value = sum(w, na.rm = TRUE), .groups = "drop")
    
    dat_pop <- pop_df %>%
      filter(year == ref_year, !is.na(simd_decile)) %>%
      mutate(
        decile_chr = as.character(simd_decile),
        source = compare_labels["population"]
      ) %>%
      group_by(source, decile_chr) %>%
      summarise(value = sum(population, na.rm = TRUE), .groups = "drop")
    
    agg <- bind_rows(dat_uw, dat_w, dat_pop) %>%
      normalise_within_source("value")
    
    if (nrow(agg) == 0) {
      warning("No rows available for SIMD population comparison: ", out_path)
      return(invisible(NULL))
    }
    
    suppressWarnings(dec_num <- as.numeric(agg$decile_chr))
    if (all(!is.na(dec_num))) {
      levs <- as.character(sort(unique(dec_num)))
    } else {
      levs <- sort(unique(agg$decile_chr))
    }
    
    y_upper <- prop_axis_upper(agg$value)
    
    agg <- agg %>%
      mutate(
        decile = factor(decile_chr, levels = levs),
        source = factor(source, levels = unname(compare_labels))
      )
    
    p <- ggplot(agg, aes(decile, value, fill = source)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, show.legend = TRUE) +
      scale_fill_manual(
        values = compare_palette,
        breaks = unname(compare_labels),
        limits = unname(compare_labels),
        drop = FALSE
      ) +
      scale_y_continuous(
        labels = percent_format(accuracy = 0.1),
        limits = c(0, y_upper),
        expand = expansion(mult = c(0, 0.02))
      ) +
      labs(x = "SIMD decile", y = "Proportion", fill = "") +
      theme_base +
      theme(axis.text.x = element_text(size = 15, colour = "black", angle = 0, hjust = 0.5))
    
    save_plot(out_path, p, w = 10, h = 5.5)
    invisible(NULL)
  }
  
  plot_simd_population_compare(
    cohorts_w$Overall,
    cohorts_ns$Overall,
    pop_clean,
    file.path(dirs$decile, "simd_2012_01_unweighted_vs_popweighted_vs_population.png")
  )
  
  # =============================================================================
  # 1b) Sex comparison:
  #     2012-01 unweighted vs 2012-01 weighted vs mid-year population 2012
  #     SHOWS PROPORTIONS
  # =============================================================================
  plot_sex_population_compare <- function(df_w_overall, df_ns_overall, pop_df, out_path) {
    if (!("sexM" %in% names(df_w_overall)) || !("sexM" %in% names(df_ns_overall))) {
      warning("sexM not found; skipping sex population comparison: ", out_path)
      return(invisible(NULL))
    }
    
    ref_month <- as.Date("2012-01-01")
    ref_year  <- 2012L
    compare_labels <- make_compare_labels(ref_month)
    compare_palette <- make_compare_palette(ref_month)
    sex_levels <- c("Females", "Males")
    
    dat_uw <- df_ns_overall %>%
      filter(date == ref_month) %>%
      mutate(
        sex_std = standardise_sex_from_cohort(sexM),
        source = compare_labels["unweighted"],
        value = 1
      ) %>%
      filter(!is.na(sex_std)) %>%
      group_by(source, sex_std) %>%
      summarise(value = sum(value), .groups = "drop")
    
    dat_w <- df_w_overall %>%
      filter(date == ref_month) %>%
      mutate(
        sex_std = standardise_sex_from_cohort(sexM),
        source = compare_labels["weighted"],
        w = get_w_from_col(., weight_col_pop)
      ) %>%
      filter(!is.na(sex_std)) %>%
      group_by(source, sex_std) %>%
      summarise(value = sum(w, na.rm = TRUE), .groups = "drop")
    
    dat_pop <- pop_df %>%
      filter(year == ref_year, !is.na(sex_std)) %>%
      mutate(source = compare_labels["population"]) %>%
      group_by(source, sex_std) %>%
      summarise(value = sum(population, na.rm = TRUE), .groups = "drop")
    
    agg <- bind_rows(dat_uw, dat_w, dat_pop) %>%
      normalise_within_source("value")
    
    if (nrow(agg) == 0) {
      warning("No rows available for sex population comparison: ", out_path)
      return(invisible(NULL))
    }
    
    y_upper <- prop_axis_upper(agg$value)
    
    agg <- agg %>%
      mutate(
        sex_std = factor(sex_std, levels = sex_levels),
        source = factor(source, levels = unname(compare_labels))
      )
    
    p <- ggplot(agg, aes(sex_std, value, fill = source)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7, show.legend = TRUE) +
      scale_fill_manual(
        values = compare_palette,
        breaks = unname(compare_labels),
        limits = unname(compare_labels),
        drop = FALSE
      ) +
      scale_y_continuous(
        labels = percent_format(accuracy = 0.1),
        limits = c(0, y_upper),
        expand = expansion(mult = c(0, 0.02))
      ) +
      labs(x = "Sex", y = "Proportion", fill = "") +
      theme_base +
      theme(axis.text.x = element_text(size = 15, colour = "black", angle = 0, hjust = 0.5))
    
    save_plot(out_path, p, w = 8, h = 5.5)
    invisible(NULL)
  }
  
  plot_sex_population_compare(
    cohorts_w$Overall,
    cohorts_ns$Overall,
    pop_clean,
    file.path(dirs$decile, "sex_2012_01_unweighted_vs_popweighted_vs_population.png")
  )
  
  # =============================================================================
  # 2) OVERVIEW PLOTS
  # =============================================================================
  make_overview <- function(out_dir_overview, cohorts, overall_prev_df, overall_label) {
    prev <- bind_rows(
      monthly_prev(overall_prev_df, overall_label),
      monthly_prev(cohorts$FE,  "FE"),
      monthly_prev(cohorts$LTC, "LTC"),
      monthly_prev(cohorts$YED, "YED")
    )
    
    prev_upper <- prop_axis_upper(prev$prop_pos)
    
    p_prev <- ggplot(prev, aes(date, prop_pos, color = cohort, group = cohort)) +
      geom_line(show.legend = TRUE) +
      geom_point(show.legend = TRUE) +
      scale_date_6m +
      scale_y_continuous(
        labels = percent_format(accuracy = 0.1),
        limits = c(0, prev_upper),
        expand = expansion(mult = c(0, 0.02))
      ) +
      scale_color_manual(
        values = cohort_palette,
        breaks = names(cohort_palette),
        limits = names(cohort_palette),
        drop = FALSE
      ) +
      labs(x = "Month", y = "Outcome Rate", color = "Cohort") +
      theme_base
    
    save_plot(file.path(out_dir_overview, "target_prevalence_by_subcohort.png"), p_prev, w = 10, h = 5.5)
    
    sizes <- bind_rows(
      monthly_size(cohorts$FE,      "FE"),
      monthly_size(cohorts$LTC,     "LTC"),
      monthly_size(cohorts$YED,     "YED"),
      monthly_size(cohorts$Overall, "Overall")
    ) %>%
      arrange(date, cohort)
    
    size_palette <- c(
      Overall = "#0072B2",
      FE      = "#D55E00",
      LTC     = "#009E73",
      YED     = "#CC79A7"
    )
    
    p_sizes <- ggplot(sizes, aes(date, n, color = cohort, group = cohort)) +
      geom_line(show.legend = TRUE) +
      geom_point(show.legend = TRUE) +
      scale_date_6m +
      scale_y_continuous(labels = comma_format()) +
      scale_color_manual(
        values = size_palette,
        breaks = names(size_palette),
        limits = names(size_palette),
        drop = FALSE
      ) +
      labs(x = "Month", y = "Patients", color = "Cohort") +
      theme_base
    
    save_plot(file.path(out_dir_overview, "cohort_sizes_by_month.png"), p_sizes, w = 10, h = 5.5)
    
    invisible(list(prev = prev, sizes = sizes))
  }
  
  overview_w <- make_overview(
    out_dir_overview = dirs$overview_w,
    cohorts          = cohorts_w,
    overall_prev_df  = cohorts_w$Overall,
    overall_label    = "Overall"
  )
  
  overview_ns <- make_overview(
    out_dir_overview = dirs$overview_ns,
    cohorts          = cohorts_ns,
    overall_prev_df  = cohorts_ns$Overall,
    overall_label    = "Overall"
  )
  
  # ---- overall unweighted vs population-weighted prevalence over time ----
  pop_prev_compare <- bind_rows(
    monthly_prev_by_weightcol(cohorts_ns$Overall, "Unweighted", wc = NULL),
    monthly_prev_by_weightcol(cohorts_w$Overall,  "Population weighted", wc = weight_col_pop)
  ) %>%
    mutate(cohort = factor(cohort, levels = c("Unweighted", "Population weighted")))
  
  pop_prev_upper <- prop_axis_upper(pop_prev_compare$prop_pos)
  
  p_pop_prev <- ggplot(pop_prev_compare, aes(date, prop_pos, color = cohort, group = cohort)) +
    geom_line(show.legend = TRUE) +
    geom_point(show.legend = TRUE) +
    scale_date_6m +
    scale_y_continuous(
      labels = percent_format(accuracy = 0.1),
      limits = c(0, pop_prev_upper),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_color_manual(
      values = c("Unweighted" = "#0072B2", "Population weighted" = "#D55E00"),
      breaks = c("Unweighted", "Population weighted"),
      limits = c("Unweighted", "Population weighted"),
      drop = FALSE
    ) +
    labs(x = "Month", y = "Outcome Rate", color = "") +
    theme_base
  
  save_plot(
    file.path(dirs$overview_pop, "target_prevalence_over_time_unweighted_vs_population_weighted.png"),
    p_pop_prev, w = 10, h = 5.5
  )
  
  # =============================================================================
  # 3) Variable overlays (weighted for age-weighted cohorts only)
  # =============================================================================
  cohorts_sub_w <- cohorts_w[c("FE","LTC","YED")]
  blacklist <- c("id","time","target","date","year_month",
                 "subcohort_FE","subcohort_LTC","subcohort_YED","subcohort_primary")
  
  vars <- setdiff(unique(unlist(lapply(cohorts_sub_w, names))), blacklist)
  
  is_num_or_log <- function(v, d) v %in% names(d) && (is.numeric(d[[v]]) || is.logical(d[[v]]))
  candidate_vars <- vars[vapply(vars, function(v) any(vapply(cohorts_sub_w, is_num_or_log, logical(1), v = v)), logical(1))]
  
  weight_vars_in_data <- vars[vapply(vars, is_analysis_weight_var, logical(1))]
  if (length(weight_vars_in_data)) {
    stale_weight_plot_files <- c(
      file.path(dirs$variables, paste0("var_", safe_filename(weight_vars_in_data), "_mean_weighted.png")),
      file.path(dirs$variables, paste0("var_", safe_filename(weight_vars_in_data), "_prev_weighted.png"))
    )
    unlink(stale_weight_plot_files[file.exists(stale_weight_plot_files)])
  }
  
  candidate_vars <- candidate_vars[!vapply(candidate_vars, is_analysis_weight_var, logical(1))]
  
  for (v in candidate_vars) {
    present <- names(cohorts_sub_w)[vapply(cohorts_sub_w, function(d) v %in% names(d), logical(1))]
    all_vals <- unlist(lapply(present, function(nm) cohorts_sub_w[[nm]][[v]]), use.names = FALSE)
    is_bin <- is_binary_vec(all_vals)
    
    dat <- bind_rows(lapply(present, function(nm) summarise_numeric(cohorts_sub_w[[nm]], v, nm)))
    y_upper_var <- if (is_bin) prop_axis_upper(dat$mean_val) else NULL
    
    p <- ggplot(dat, aes(date, mean_val, color = cohort, group = cohort)) +
      geom_line(show.legend = TRUE) +
      geom_point(show.legend = TRUE) +
      scale_date_6m +
      scale_color_manual(
        values = sub_palette,
        breaks = names(sub_palette),
        limits = names(sub_palette),
        drop = FALSE
      ) +
      labs(x = "Month", y = if (is_bin) "Weighted prevalence" else "Weighted mean", color = "Cohort") +
      theme_base +
      (if (is_bin)
        scale_y_continuous(
          labels = percent_format(accuracy = 0.1),
          limits = c(0, y_upper_var),
          expand = expansion(mult = c(0, 0.02))
        )
       else
         scale_y_continuous(labels = label_number(accuracy = 0.01))
      )
    
    suffix <- if (is_bin) "_prev_weighted.png" else "_mean_weighted.png"
    save_plot(file.path(dirs$variables, paste0("var_", safe_filename(v), suffix)), p, w = 10, h = 5.5)
  }
  
  # =============================================================================
  # 4) Age overlays (Jan-2012 vs Oct-2022) for age-weighted + unweighted
  # =============================================================================
  date1 <- as.Date("2012-01-01")
  date2 <- as.Date("2022-10-01")
  period_palette <- c(
    "2012-01" = "#0072B2",
    "2022-10" = "#D55E00"
  )
  
  plot_age_overlay <- function(df, cohort_name, out_path, weighted = FALSE) {
    if (!("age" %in% names(df))) return(invisible(NULL))
    
    dat <- df %>%
      filter(date %in% c(date1, date2), !is.na(age)) %>%
      mutate(
        period = factor(ifelse(date == date1, "2012-01", "2022-10"),
                        levels = c("2012-01","2022-10")),
        w = if (weighted) get_w(., cohort = cohort_name) else 1
      )
    
    if (nrow(dat) == 0) return(invisible(NULL))
    
    p <- ggplot(dat, aes(x = age, fill = period, weight = w)) +
      geom_density(alpha = 0.35, show.legend = TRUE) +
      scale_fill_manual(
        values = period_palette,
        breaks = c("2012-01", "2022-10"),
        limits = c("2012-01", "2022-10"),
        drop = FALSE
      ) +
      scale_x_continuous(limits = c(0, NA)) +
      labs(x = "Age", y = "Density", fill = "") +
      theme_base
    
    save_plot(out_path, p, w = 9, h = 5)
  }
  
  for (nm in names(cohorts_w)) {
    plot_age_overlay(
      cohorts_w[[nm]],
      cohort_name = nm,
      out_path = file.path(dirs$age_w, paste0("age_overlay_", nm, "_weighted.png")),
      weighted = TRUE
    )
  }
  
  for (nm in names(cohorts_ns)) {
    plot_age_overlay(
      cohorts_ns[[nm]],
      cohort_name = nm,
      out_path = file.path(dirs$age_ns, paste0("age_overlay_", nm, "_unweighted.png")),
      weighted = FALSE
    )
  }
  
  # =============================================================================
  # 5) Population-weighted age density comparison:
  #    2012-01 unweighted vs 2012-01 weighted vs mid-year population 2012
  #    2022-10 unweighted vs 2022-10 weighted vs mid-year population 2022
  # =============================================================================
  plot_age_population_compare <- function(df_w_overall, df_ns_overall, pop_df, ref_month, out_path,
                                          min_age_included = MIN_AGE_INCLUDED) {
    ref_year <- year(ref_month)
    compare_labels <- make_compare_labels(ref_month)
    compare_palette <- make_compare_palette(ref_month)
    
    dat_uw <- df_ns_overall %>%
      filter(date == ref_month, !is.na(age)) %>%
      transmute(
        age = suppressWarnings(as.numeric(age)),
        source = compare_labels["unweighted"],
        w = 1
      ) %>%
      filter(!is.na(age), age >= min_age_included)
    
    dat_w_src <- df_w_overall %>%
      filter(date == ref_month, !is.na(age))
    
    dat_w <- dat_w_src %>%
      transmute(
        age = suppressWarnings(as.numeric(age)),
        source = compare_labels["weighted"],
        w = get_w_from_col(dat_w_src, weight_col_pop)
      ) %>%
      filter(!is.na(age), age >= min_age_included)
    
    dat_pop <- pop_df %>%
      filter(year == ref_year, !is.na(age_num), age_num >= min_age_included) %>%
      transmute(
        age = age_num,
        source = compare_labels["population"],
        w = population
      )
    
    dat <- bind_rows(dat_uw, dat_w, dat_pop) %>%
      mutate(source = factor(source, levels = unname(compare_labels)))
    
    if (nrow(dat) == 0) {
      warning("No rows available for population age-density comparison: ", out_path)
      return(invisible(NULL))
    }
    
    p <- ggplot(dat, aes(x = age, color = source, fill = source, weight = w)) +
      geom_density(alpha = 0.20, show.legend = TRUE) +
      scale_color_manual(
        values = compare_palette,
        breaks = unname(compare_labels),
        limits = unname(compare_labels),
        drop = FALSE
      ) +
      scale_fill_manual(
        values = compare_palette,
        breaks = unname(compare_labels),
        limits = unname(compare_labels),
        drop = FALSE
      ) +
      scale_x_continuous(limits = c(min_age_included, NA)) +
      labs(x = "Age", y = "Density", color = "", fill = "") +
      theme_base +
      theme(axis.text.x = element_text(size = 15, colour = "black", angle = 0, hjust = 0.5))
    
    save_plot(out_path, p, w = 10, h = 5.5)
    invisible(NULL)
  }
  
  plot_age_population_compare(
    cohorts_w$Overall,
    cohorts_ns$Overall,
    pop_clean,
    ref_month = as.Date("2012-01-01"),
    out_path = file.path(dirs$age_pop, "age_density_2012_01_unweighted_vs_popweighted_vs_population.png")
  )
  
  plot_age_population_compare(
    cohorts_w$Overall,
    cohorts_ns$Overall,
    pop_clean,
    ref_month = as.Date("2022-10-01"),
    out_path = file.path(dirs$age_pop, "age_density_2022_10_unweighted_vs_popweighted_vs_population.png")
  )
  
  # =============================================================================
  # 6) Target reasons (quarterly stacked) — grouped only
  # =============================================================================
  if (!"first_event_reason" %in% names(df_all)) {
    warning("first_event_reason not found; skipping target reason plots.")
    return(invisible(list(
      FE = cohorts_w$FE,
      LTC = cohorts_w$LTC,
      YED = cohorts_w$YED,
      Overall = cohorts_w$Overall,
      prev = overview_w$prev,
      sizes = overview_w$sizes,
      pop_prev_compare = pop_prev_compare,
      weight_diagnostics = weight_diagnostics
    )))
  }
  
  prep_quarter_grouped <- function(df, cohort, weighted = FALSE, wc = NULL, min_quarter_count = 5) {
    w <- if (weighted) get_w_from_col(df, wc %||% resolve_weight_col(cohort)) else rep(1, nrow(df))
    other_reason <- "OTHER_OR_UNMAPPED"
    valid_levels <- grouped_levels
    quarter_dat <- df %>%
      mutate(w = w) %>%
      filter(!is.na(target), target == 1L) %>%
      mutate(
        quarter = floor_date(date, unit = "quarter"),
        reason = ifelse(
          is.na(first_event_reason) | first_event_reason == "",
          other_reason,
          first_event_reason
        ),
        reason = ifelse(reason %in% valid_levels, reason, other_reason)
      ) %>%
      group_by(quarter, reason) %>%
      summarise(value = sum(w, na.rm = TRUE), .groups = "drop")
    
    if (nrow(quarter_dat) == 0) {
      return(tibble(
        quarter = as.Date(character()),
        reason = factor(character(), levels = grouped_levels),
        value = numeric(),
        total = numeric(),
        prop = numeric(),
        cohort = character()
      ))
    }
    
    all_quarters <- sort(unique(quarter_dat$quarter))
    if (length(all_quarters) > 1) {
      all_quarters <- all_quarters[all_quarters < max(all_quarters, na.rm = TRUE)]
    } else {
      all_quarters <- as.Date(character())
    }
    
    if (length(all_quarters) == 0) {
      return(tibble(
        quarter = as.Date(character()),
        reason = factor(character(), levels = grouped_levels),
        value = numeric(),
        total = numeric(),
        prop = numeric(),
        cohort = character()
      ))
    }
    
    quarter_dat <- quarter_dat %>%
      filter(quarter %in% all_quarters)
    
    expanded <- tidyr::expand_grid(
      quarter = all_quarters,
      reason = valid_levels
    ) %>%
      left_join(quarter_dat, by = c("quarter", "reason")) %>%
      mutate(value = replace_na(value, 0))
    
    rare_reasons <- expanded %>%
      filter(reason != other_reason) %>%
      group_by(reason) %>%
      summarise(any_lt5 = any(value < min_quarter_count), .groups = "drop") %>%
      filter(any_lt5) %>%
      pull(reason)
    
    collapsed <- expanded %>%
      mutate(reason = ifelse(reason %in% rare_reasons, other_reason, reason)) %>%
      group_by(quarter, reason) %>%
      summarise(value = sum(value, na.rm = TRUE), .groups = "drop")
    
    other_bad <- collapsed %>%
      filter(reason == other_reason) %>%
      summarise(any_lt5 = any(value < min_quarter_count), .groups = "drop") %>%
      pull(any_lt5)
    
    if (length(other_bad) == 0) other_bad <- FALSE
    
    if (isTRUE(other_bad)) {
      collapsed <- collapsed %>%
        filter(reason != other_reason)
    }
    
    keep_levels <- grouped_levels[grouped_levels %in% unique(collapsed$reason)]
    
    collapsed %>%
      mutate(reason = factor(reason, levels = keep_levels)) %>%
      group_by(quarter) %>%
      mutate(total = sum(value), prop = ifelse(total > 0, value / total, NA_real_)) %>%
      ungroup() %>%
      mutate(cohort = cohort)
  }
  
  plot_stack <- function(dat, y, palette, levels, ylab) {
    levels_use <- levels[levels %in% unique(as.character(dat$reason))]
    p <- ggplot(dat, aes(quarter, .data[[y]], fill = reason)) +
      geom_col(show.legend = TRUE) +
      scale_fill_manual(
        values = palette[levels_use],
        breaks = levels_use,
        limits = levels_use,
        labels = grouped_labels[levels_use],
        drop = TRUE,
        na.translate = FALSE
      ) +
      scale_date_6m +
      labs(x = "Quarter", y = ylab, fill = "ICD-10 chapter") +
      theme_base +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)
      )
    
    if (y == "prop") {
      p <- p + scale_y_continuous(
        labels = percent_format(accuracy = 0.1),
        limits = c(0, 1),
        expand = expansion(mult = c(0, 0.02))
      )
    }
    if (y == "value") {
      p <- p + scale_y_continuous(labels = comma_format())
    }
    p
  }
  
  plot_reason_lines <- function(dat, y, palette, levels, ylab) {
    levels_use <- levels[levels %in% unique(as.character(dat$reason))]
    p <- ggplot(
      dat,
      aes(
        x = quarter,
        y = .data[[y]],
        color = reason,
        group = reason
      )
    ) +
      geom_line(linewidth = 0.8, show.legend = TRUE, na.rm = TRUE) +
      geom_point(size = 1.6, show.legend = FALSE, na.rm = TRUE) +
      scale_color_manual(
        values = palette[levels_use],
        breaks = levels_use,
        limits = levels_use,
        labels = grouped_labels[levels_use],
        drop = TRUE,
        na.translate = FALSE
      ) +
      scale_date_6m +
      labs(x = "Quarter", y = ylab, color = "ICD-10 chapter") +
      theme_base +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)
      )
    
    if (y == "prop") {
      p <- p + scale_y_continuous(
        labels = percent_format(accuracy = 0.1),
        limits = c(0, 1),
        expand = expansion(mult = c(0, 0.02))
      )
    }
    if (y == "value") {
      p <- p + scale_y_continuous(labels = comma_format())
    }
    p
  }
  
  emit_grouped_reason_plots <- function(cohort_list, out_group, suffix = "", weighted = FALSE) {
    for (nm in names(cohort_list)) {
      dfc <- cohort_list[[nm]]
      if (is.null(dfc) || nrow(dfc) == 0) next
      
      qg <- prep_quarter_grouped(dfc, nm, weighted = weighted)
      if (nrow(qg) == 0) next
      count_label <- if (weighted) "Weighted count" else "Unweighted count"
      
      save_plot(
        file.path(out_group, paste0("targets_reason_", nm, "_prop", suffix, ".png")),
        plot_stack(qg, "prop", grouped_palette, grouped_levels, "Proportion"),
        w = 11, h = 6
      )
      save_plot(
        file.path(out_group, paste0("targets_reason_", nm, "_counts", suffix, ".png")),
        plot_stack(qg, "value", grouped_palette, grouped_levels, count_label),
        w = 11, h = 6
      )
      save_plot(
        file.path(out_group, paste0("targets_reason_", nm, "_counts_line", suffix, ".png")),
        plot_reason_lines(qg, "value", grouped_palette, grouped_levels, count_label),
        w = 11, h = 6
      )
    }
  }
  
  emit_grouped_reason_plots(
    cohorts_w,
    dirs$tr_w_g,
    suffix = "_weighted",
    weighted = TRUE
  )
  
  emit_grouped_reason_plots(
    cohorts_ns,
    dirs$tr_ns_g,
    suffix = "_unweighted",
    weighted = FALSE
  )
  
  # ---- grouped target reasons for population-weighted OVERALL only ----
  qg_pop_overall <- prep_quarter_grouped(
    cohorts_w$Overall,
    cohort = "Overall",
    weighted = TRUE,
    wc = weight_col_pop
  )
  
  save_plot(
    file.path(dirs$tr_pop_g, "targets_reason_Overall_prop_population_weighted.png"),
    plot_stack(qg_pop_overall, "prop", grouped_palette, grouped_levels, "Proportion"),
    w = 11, h = 6
  )
  
  save_plot(
    file.path(dirs$tr_pop_g, "targets_reason_Overall_counts_population_weighted.png"),
    plot_stack(qg_pop_overall, "value", grouped_palette, grouped_levels, "Weighted count"),
    w = 11, h = 6
  )
  
  save_plot(
    file.path(dirs$tr_pop_g, "targets_reason_Overall_counts_line_population_weighted.png"),
    plot_reason_lines(qg_pop_overall, "value", grouped_palette, grouped_levels, "Weighted count"),
    w = 11, h = 6
  )
  
  invisible(list(
    FE = cohorts_w$FE,
    LTC = cohorts_w$LTC,
    YED = cohorts_w$YED,
    Overall = cohorts_w$Overall,
    prev = overview_w$prev,
    sizes = overview_w$sizes,
    pop_prev_compare = pop_prev_compare,
    weight_diagnostics = weight_diagnostics
  ))
}

`%||%` <- function(x, y) if (is.null(x)) y else x