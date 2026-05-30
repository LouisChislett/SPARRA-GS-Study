# =============================================================================
# weightedDataDrift.R
#
# Author      : Louis Chislett
#
# Purpose     :
#   Weighted descriptive data-drift analysis for SPARRA CPM v3 RAW cohorts.
#
#   Key design choices:
#     - Produces DESCRIPTIVE distances only:
#         * weighted KS distance for continuous variables
#         * weighted total variation (TV) distance for categorical variables
# =============================================================================

weightedDataDrift <- function(
    df_all = NULL,
    in_file = NULL,
    month1,
    month2,
    out_dir = "data drift",
    include_overall = FALSE,
    seed = 123,
    make_plots = TRUE,
    plot_top_n = 10,
    continuous_vars = NULL,
    categorical_vars = NULL,
    write_outputs = TRUE,
    summary_csv_name = "resultsdataDrift.csv",
    categorical_props_csv_name = "weighted_drift_categorical_proportions.csv",
    cohort_weight_map = c(
      FE  = "w_age_FE_trim",
      LTC = "w_age_LTC_trim",
      YED = "w_age_YED_trim",
      Overall = "w_age_trim"
    )
) {
  # Load dependencies locally inside the function so the script can be sourced
  # without requiring the caller to attach these packages first. Startup
  # messages are suppressed to keep batch logs cleaner.
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(purrr)
    library(lubridate)
    library(ggplot2)
    library(scales)
    library(stringr)
    library(forcats)
    library(fst)
  })
  
  # Set a seed for reproducibility. The current calculations are deterministic,
  # but keeping this here makes future sampling/plot extensions reproducible.
  set.seed(seed)
  
  # Small null-coalescing helper: return y when x is NULL or empty.
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
  
  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  # Create a directory if it does not already exist.
  ensure_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
  }
  
  # Convert variable names into filesystem-safe names for plot files.
  safe_filename <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)
  
  # Try to determine where this script lives. This supports relative input paths
  # both when the script is run with Rscript and when it is sourced.
  get_script_dir <- function() {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cmd_args, value = TRUE)
    if (length(file_arg)) {
      return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
    }
    
    ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
    }
    
    getwd()
  }
  
  script_dir <- get_script_dir()
  
  # Resolve a relative path against the current working directory and the script
  # directory. If no matching file is found, return the original path so the
  # downstream error message remains informative.
  resolve_local_path <- function(path) {
    if (is.null(path) || !length(path) || is.na(path) || !nzchar(path)) return(path)
    if (file.exists(path)) return(path)
    
    candidates <- unique(c(
      file.path(getwd(), path),
      file.path(script_dir, path)
    ))
    hit <- candidates[file.exists(candidates)][1]
    if (is.na(hit) || !length(hit)) return(path)
    hit
  }
  
  # Standardise binary target encodings to integer 0/1. This handles logical,
  # numeric, factor, and common character encodings such as TRUE/FALSE and YES/NO.
  to01 <- function(x) {
    if (is.logical(x)) return(as.integer(x))
    
    if (is.numeric(x)) {
      return(as.integer(ifelse(is.na(x), NA, round(x))))
    }
    
    if (is.factor(x)) x <- as.character(x)
    
    if (is.character(x)) {
      tx <- toupper(trimws(x))
      out <- ifelse(tx %in% c("1", "TRUE", "T", "YES", "Y"), 1L,
                    ifelse(tx %in% c("0", "FALSE", "F", "NO", "N"), 0L, NA_integer_))
      return(out)
    }
    
    stop("Unsupported target type: ", paste(class(x), collapse = ", "))
  }
  
  # Convert month inputs to the first day of the relevant month. This gives a
  # stable Date value for exact matching against the generated `date` column.
  to_date_month <- function(x) {
    if (inherits(x, "Date")) return(as.Date(format(x, "%Y-%m-01")))
    if (inherits(x, "POSIXt")) return(as.Date(format(as.Date(x), "%Y-%m-01")))
    x <- as.character(x)
    if (grepl("^\\d{4}-\\d{2}$", x[1])) return(as.Date(paste0(x, "-01")))
    as.Date(x)
  }
  
  # Ensure every cohort has both `year_month` and `date` columns, deriving them
  # from `time` if necessary, then sort chronologically.
  add_dates <- function(df) {
    if (!nrow(df)) return(df)
    df %>%
      mutate(
        year_month = if (!"year_month" %in% names(.)) format(as_datetime(time, tz = "UTC"), "%Y-%m") else year_month,
        date       = if (!"date" %in% names(.)) as.Date(paste0(year_month, "-01")) else as.Date(date)
      ) %>%
      arrange(date)
  }
  
  # Coerce weights to numeric and set unusable values to zero. Negative, missing,
  # infinite, and non-numeric weights should not contribute to weighted summaries.
  sanitize_weights <- function(w) {
    w <- suppressWarnings(as.numeric(w))
    w[!is.finite(w) | is.na(w) | w < 0] <- 0
    w
  }
  
  # Pull a cohort-specific weight column if present; otherwise fall back to
  # equal weighting.
  get_w_from_col <- function(df, wc) {
    if (!is.null(wc) && wc %in% names(df)) {
      return(sanitize_weights(df[[wc]]))
    }
    rep(1, nrow(df))
  }
  
  # Compute weighted quantiles by sorting values and finding the first point
  # where cumulative weight reaches each requested probability.
  weighted_quantile <- function(x, w, probs = c(0.25, 0.5, 0.75)) {
    ok <- is.finite(x) & !is.na(x) & is.finite(w) & !is.na(w) & w > 0
    x <- as.numeric(x[ok])
    w <- as.numeric(w[ok])
    if (!length(x)) return(rep(NA_real_, length(probs)))
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    cw <- cumsum(w) / sum(w)
    vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
  }
  
  # Identify numeric/logical variables that should be treated as categorical
  # because they only contain binary 0/1-style values.
  is_binary_like <- function(x) {
    ux <- sort(unique(x[!is.na(x)]))
    length(ux) > 0 && length(ux) <= 2 && all(ux %in% c(0, 1, FALSE, TRUE))
  }
  
  # Infer whether a variable should be analysed as continuous or categorical.
  # Character, factor, logical, and binary-like numeric variables are categorical;
  # other numeric variables are continuous.
  infer_variable_type <- function(x) {
    x_no_na <- x[!is.na(x)]
    if (!length(x_no_na)) return("unknown")
    
    if (is.factor(x) || is.character(x)) return("categorical")
    if (is.logical(x)) return("categorical")
    
    if (is.numeric(x)) {
      if (is_binary_like(x)) return("categorical")
      return("continuous")
    }
    
    "categorical"
  }
  
  # Build a weighted empirical CDF table for a continuous variable. Duplicate
  # values are collapsed by summing their weights before calculating the CDF.
  weighted_ecdf_table <- function(x, w) {
    ok <- is.finite(x) & !is.na(x) & is.finite(w) & !is.na(w) & w > 0
    x <- as.numeric(x[ok])
    w <- as.numeric(w[ok])
    if (!length(x)) return(tibble(value = numeric(), cdf = numeric()))
    tibble(value = x, w = w) %>%
      group_by(value) %>%
      summarise(w = sum(w), .groups = "drop") %>%
      arrange(value) %>%
      mutate(cdf = cumsum(w) / sum(w)) %>%
      select(value, cdf)
  }
  
  # Weighted Kolmogorov-Smirnov distance for continuous variables. This is the
  # largest absolute gap between the two weighted empirical CDFs.
  weighted_ks_distance <- function(x1, w1, x2, w2) {
    t1 <- weighted_ecdf_table(x1, w1)
    t2 <- weighted_ecdf_table(x2, w2)
    if (!nrow(t1) || !nrow(t2)) return(NA_real_)
    
    support <- sort(unique(c(t1$value, t2$value)))
    
    step_eval <- function(tbl, s) {
      idx <- findInterval(s, tbl$value)
      out <- numeric(length(s))
      ok <- idx > 0
      out[ok] <- tbl$cdf[idx[ok]]
      out
    }
    
    f1 <- step_eval(t1, support)
    f2 <- step_eval(t2, support)
    max(abs(f1 - f2), na.rm = TRUE)
  }
  
  # Weighted total variation distance for categorical variables. This compares
  # weighted category proportions across the union of observed levels.
  weighted_tv_distance <- function(x1, w1, x2, w2) {
    ok1 <- !is.na(x1) & !is.na(w1) & is.finite(w1) & w1 > 0
    ok2 <- !is.na(x2) & !is.na(w2) & is.finite(w2) & w2 > 0
    
    x1 <- as.character(x1[ok1])
    x2 <- as.character(x2[ok2])
    w1 <- as.numeric(w1[ok1])
    w2 <- as.numeric(w2[ok2])
    
    if (!length(x1) || !length(x2)) return(NA_real_)
    
    levs <- sort(unique(c(x1, x2)))
    
    p1 <- tibble(level = x1, w = w1) %>%
      group_by(level) %>%
      summarise(p = sum(w), .groups = "drop") %>%
      mutate(p = p / sum(p))
    
    p2 <- tibble(level = x2, w = w2) %>%
      group_by(level) %>%
      summarise(p = sum(w), .groups = "drop") %>%
      mutate(p = p / sum(p))
    
    full <- tibble(level = levs) %>%
      left_join(p1, by = "level") %>%
      rename(p1 = p) %>%
      left_join(p2, by = "level") %>%
      rename(p2 = p) %>%
      mutate(
        p1 = coalesce(p1, 0),
        p2 = coalesce(p2, 0)
      )
    
    0.5 * sum(abs(full$p1 - full$p2))
  }
  
  # Produce weighted descriptive statistics for a continuous variable in one
  # month/cohort: count, weight sum, mean, quartiles, and range.
  summarise_continuous <- function(x, w) {
    ok <- is.finite(x) & !is.na(x) & is.finite(w) & !is.na(w) & w > 0
    x <- as.numeric(x[ok])
    w <- as.numeric(w[ok])
    if (!length(x)) {
      return(tibble(
        n_non_missing = 0L,
        weight_sum = 0,
        mean_w = NA_real_,
        q25_w = NA_real_,
        median_w = NA_real_,
        q75_w = NA_real_,
        min = NA_real_,
        max = NA_real_
      ))
    }
    qs <- weighted_quantile(x, w, c(0.25, 0.5, 0.75))
    tibble(
      n_non_missing = length(x),
      weight_sum = sum(w),
      mean_w = sum(x * w) / sum(w),
      q25_w = qs[1],
      median_w = qs[2],
      q75_w = qs[3],
      min = min(x),
      max = max(x)
    )
  }
  
  # Produce weighted descriptive statistics for a categorical variable in one
  # month/cohort: count, weight sum, number of levels, top level, and level
  # proportions.
  summarise_categorical <- function(x, w) {
    ok <- !is.na(x) & !is.na(w) & is.finite(w) & w > 0
    x <- as.character(x[ok])
    w <- as.numeric(w[ok])
    if (!length(x)) {
      return(list(
        summary = tibble(n_non_missing = 0L, weight_sum = 0, n_levels = 0L, top_level = NA_character_, top_prop_w = NA_real_),
        proportions = tibble(level = character(), prop_w = numeric())
      ))
    }
    
    props <- tibble(level = x, w = w) %>%
      group_by(level) %>%
      summarise(weight = sum(w), .groups = "drop") %>%
      mutate(prop_w = weight / sum(weight)) %>%
      arrange(desc(prop_w), level)
    
    list(
      summary = tibble(
        n_non_missing = length(x),
        weight_sum = sum(w),
        n_levels = nrow(props),
        top_level = props$level[1],
        top_prop_w = props$prop_w[1]
      ),
      proportions = props %>% select(level, prop_w)
    )
  }
  
  # ---- plotting helpers ----
  palette_months <- c(
    reference = "#0072B2",  # blue
    comparison = "#D55E00"  # vermilion
  )
  
  base_theme <- theme_minimal(base_size = 20) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(size = 16),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 16),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
  
  # Save an overlaid weighted histogram/bar chart for a continuous variable.
  # Integer-like variables are plotted as discrete bars; other variables are
  # binned into equal-width intervals.
  plot_continuous_histogram <- function(df1, df2, var, month1_lab, month2_lab, out_path, histogram_bins = 30) {
    d1 <- df1 %>%
      transmute(value = .data[[var]], w = .data[["..w.."]], source = month1_lab)
    d2 <- df2 %>%
      transmute(value = .data[[var]], w = .data[["..w.."]], source = month2_lab)
    
    hist_dat <- bind_rows(d1, d2) %>%
      filter(is.finite(value), !is.na(value), is.finite(w), !is.na(w), w > 0)
    if (!nrow(hist_dat)) return(invisible(NULL))
    
    xlab <- str_replace_all(var, "_", " ")
    is_integer_like <- all(abs(hist_dat$value - round(hist_dat$value)) < 1e-8)
    
    if (is_integer_like) {
      plot_dat <- hist_dat %>%
        mutate(bin = as.integer(round(value))) %>%
        group_by(source, bin) %>%
        summarise(weight = sum(w), .groups = "drop") %>%
        group_by(source) %>%
        mutate(prop = weight / sum(weight)) %>%
        ungroup()
      
      p <- ggplot(plot_dat, aes(x = factor(bin), y = prop, fill = source)) +
        geom_col(
          position = "identity",
          alpha = 0.45,
          colour = "white",
          linewidth = 0.2,
          width = 0.9
        )
    } else {
      x_min <- min(hist_dat$value, na.rm = TRUE)
      x_max <- max(hist_dat$value, na.rm = TRUE)
      
      if (!is.finite(x_min) || !is.finite(x_max)) return(invisible(NULL))
      
      if (x_min == x_max) {
        breaks <- c(x_min - 0.5, x_max + 0.5)
      } else {
        breaks <- seq(x_min, x_max, length.out = histogram_bins + 1)
      }
      
      bin_template <- cut(hist_dat$value, breaks = breaks, include.lowest = TRUE, right = FALSE)
      bin_levels <- levels(bin_template)
      mids <- head(breaks, -1) + diff(breaks) / 2
      widths <- diff(breaks)
      
      plot_dat <- hist_dat %>%
        mutate(bin = cut(value, breaks = breaks, include.lowest = TRUE, right = FALSE)) %>%
        group_by(source, bin) %>%
        summarise(weight = sum(w), .groups = "drop") %>%
        tidyr::complete(
          source = c(month1_lab, month2_lab),
          bin = factor(bin_levels, levels = bin_levels),
          fill = list(weight = 0)
        ) %>%
        group_by(source) %>%
        mutate(prop = weight / sum(weight)) %>%
        ungroup() %>%
        mutate(
          mid = mids[match(as.character(bin), bin_levels)],
          width = widths[match(as.character(bin), bin_levels)]
        )
      
      p <- ggplot(plot_dat, aes(x = mid, y = prop, fill = source)) +
        geom_col(
          aes(width = width),
          position = "identity",
          alpha = 0.45,
          colour = "white",
          linewidth = 0.2
        )
    }
    
    p <- p +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      scale_fill_manual(
        values = setNames(unname(palette_months), c(month1_lab, month2_lab)),
        breaks = c(month1_lab, month2_lab)
      ) +
      labs(
        x = xlab,
        y = "Weighted proportion",
        fill = "Month"
      ) +
      base_theme
    
    ggsave(out_path, p, width = 9, height = 5.5, dpi = 300, bg = "white")
    invisible(NULL)
  }
  
  
  # Save side-by-side weighted category proportions for a categorical variable.
  # `distance_value` is accepted for interface symmetry, although the current plot
  # does not print it in the figure.
  plot_categorical_compare <- function(df1, df2, var, month1_lab, month2_lab, distance_value, out_path) {
    s1 <- summarise_categorical(df1[[var]], df1[["..w.."]])$proportions %>% mutate(source = month1_lab)
    s2 <- summarise_categorical(df2[[var]], df2[["..w.."]])$proportions %>% mutate(source = month2_lab)
    pdat <- bind_rows(s1, s2)
    if (!nrow(pdat)) return(invisible(NULL))
    
    pdat <- pdat %>%
      mutate(level = fct_reorder(level, prop_w, .fun = max, .desc = TRUE))
    
    xlab <- str_replace_all(var, "_", " ")
    
    p <- ggplot(pdat, aes(x = level, y = prop_w, fill = source)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      scale_y_continuous(labels = percent_format(accuracy = 1)) +
      scale_fill_manual(
        values = setNames(unname(palette_months), c(month1_lab, month2_lab)),
        breaks = c(month1_lab, month2_lab)
      ) +
      labs(
        y = "Weighted proportion",
        fill = "Month"
      ) +
      base_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16))
    
    ggsave(out_path, p, width = 10, height = 6, dpi = 300, bg = "white")
    invisible(NULL)
  }
  
  # Canonical cohort-specific model inputs, aligned to models.R::model_def().
  # These are the default RAW variables analysed for each cohort unless the user
  # supplies `continuous_vars` and/or `categorical_vars` to restrict the analysis.
  default_model_vars_by_cohort <- list(
    FE = c(
      "num_emergency_admissions","emergency_bed_days","num_alcohol_admissions",
      "num_ae2_attendances","num_elective_admissions","num_outpatient_appointment_general",
      "num_bnf_sections","pis_incontinence","pis_respiratory","pis_cns","pis_infections",
      "pis_endocrine","parkinsons_indicated","numLTCs_resulting_in_admin","age","decile"
    ),
    LTC = c(
      "num_emergency_admissions","emergency_bed_days","num_alcohol_drug_admissions",
      "num_ae2_attendances","num_daycase_admissions","num_elective_admissions",
      "elective_bed_days","num_psych_admissions","num_outpatient_appointment_general",
      "numLTCs_resulting_in_admin","num_bnf_sections","epilepsy_indicated","MS_indicated",
      "parkinsons_indicated","pis_infections","pis_sub_depend","pis_dementia",
      "pis_corticosteroids","pis_fluids","pis_nutrition","pis_vitamins","pis_bandages",
      "pis_catheters","pis_stoma","decile"
    ),
    YED = c(
      "num_alcohol_drug_admissions","num_emergency_admissions","emergency_bed_days",
      "num_ae2_attendances","num_psych_admissions","num_elective_admissions",
      "num_outpatient_appointment_general","num_outpatient_appointment_psych",
      "numLTCs_resulting_in_admin","num_bnf_sections","pis_gut_motility","pis_antisecretory",
      "pis_intestinal","pis_antifibrinolytic","pis_anticoagulant","pis_stoma","pis_mucolytics",
      "pis_diabetes","pis_corticosteroids","pis_fluids","pis_vitamins","pis_cns","dementia_indicated"
    )
  )
  
  # ---------------------------------------------------------------------------
  # Input loading
  # Either use the supplied data frame or read an fst file from disk.
  # ---------------------------------------------------------------------------
  if (is.null(df_all)) {
    if (is.null(in_file) || !length(in_file) || is.na(in_file) || !nzchar(in_file)) {
      stop("weightedDataDrift() requires either df_all or a valid in_file.")
    }
    
    # Support relative paths by checking both the working directory and the
    # directory containing this script.
    in_file_resolved <- resolve_local_path(in_file)
    if (!file.exists(in_file_resolved)) {
      stop(sprintf("weightedDataDrift() could not find input file: %s", in_file))
    }
    
    message("[weightedDataDrift] Reading fst input: ", in_file_resolved)
    df_all <- read_fst(in_file_resolved)
  }
  
  # Convert the target to a consistent 0/1 representation before cohort
  # filtering. If the target is absent, continue with a warning because the
  # cohort builder may still fail later depending on its requirements.
  if ("target" %in% names(df_all)) {
    df_all <- df_all %>% mutate(target = to01(target))
  } else {
    warning("No 'target' column found in df_all; drift cohort filtering may fail.")
  }
  
  df_all <- add_dates(df_all)
  
  # ---------------------------------------------------------------------------
  # Build weighted RAW cohorts
  # This function deliberately uses the RAW cohort builder so drift is assessed
  # before CPM binning/transformation.
  # ---------------------------------------------------------------------------
  if (!exists("make_sparra_v3_datasets_raw", mode = "function")) {
    stop("make_sparra_v3_datasets_raw() not found. Source cohortBuilder.R first.")
  }
  
  # Normalise comparison months and create YYYY-MM labels for outputs.
  month1 <- to_date_month(month1)
  month2 <- to_date_month(month2)
  month1_lab <- format(month1, "%Y-%m")
  month2_lab <- format(month2, "%Y-%m")
  
  # Prepare the output folder structure up front, including separate folders for
  # continuous histograms and categorical plots.
  ensure_dir(out_dir)
  ensure_dir(file.path(out_dir, "plots"))
  ensure_dir(file.path(out_dir, "plots", "continuous"))
  ensure_dir(file.path(out_dir, "plots", "continuous", "histogram"))
  ensure_dir(file.path(out_dir, "plots", "categorical"))
  
  # Build the FE, LTC, and YED raw cohorts, then keep only rows with a known
  # target and ensure each cohort has date/month columns.
  cohorts <- make_sparra_v3_datasets_raw(df_all)
  cohorts <- lapply(cohorts, function(x) add_dates(x %>% filter(!is.na(target))))
  
  # Optionally add an Overall cohort by row-binding the three cohort datasets
  # and removing duplicate rows.
  if (isTRUE(include_overall)) {
    keep_cols <- unique(unlist(lapply(cohorts, names)))
    cohorts$Overall <- bind_rows(
      cohorts$FE  %>% select(any_of(keep_cols)),
      cohorts$LTC %>% select(any_of(keep_cols)),
      cohorts$YED %>% select(any_of(keep_cols))
    ) %>%
      distinct() %>%
      arrange(date)
  }
  
  # Only analyse cohorts that also have an entry in the weight map.
  cohort_names <- names(cohorts)
  cohort_names <- cohort_names[cohort_names %in% names(cohort_weight_map)]
  
  # ---------------------------------------------------------------------------
  # Variable specification
  # Decide which variables are eligible for analysis and infer how each should be
  # handled. Operational columns and weight columns are excluded from drift tests.
  # ---------------------------------------------------------------------------
  # Collect all variables present across selected cohorts. Kept mainly as a
  # useful diagnostic object during debugging.
  all_present_vars <- sort(unique(unlist(lapply(cohorts[cohort_names], names))))
  # Columns that should never be assessed as predictors for drift. This removes
  # identifiers, dates, targets, cohort flags, and weights.
  blacklist <- c(
    "id","time","target","date","year_month",
    "subcohort_FE","subcohort_LTC","subcohort_YED","subcohort_primary",
    "age_band_std",
    "w_age_raw","w_age_stab","w_age_trim",
    "w_age_FE_raw","w_age_FE_stab","w_age_FE_trim",
    "w_age_LTC_raw","w_age_LTC_stab","w_age_LTC_trim",
    "w_age_YED_raw","w_age_YED_stab","w_age_YED_trim",
    "w_pop_raw","w_pop_stab","w_pop_trim"
  )
  
  `%notin%` <- Negate(`%in%`)
  
  # Infer the analysis type for the requested variables in a particular cohort.
  infer_types_for_vars <- function(df, vars) {
    vars <- intersect(unique(vars), names(df))
    if (!length(vars)) {
      return(tibble(variable = character(), variable_type = character()))
    }
    
    tibble(
      variable = vars,
      variable_type = purrr::map_chr(vars, function(v) infer_variable_type(df[[v]]))
    )
  }
  
  # Start from the canonical model variable list for each cohort, remove
  # blacklisted columns, and keep only variables that actually exist in that
  # cohort dataset.
  model_vars_by_cohort <- lapply(cohort_names, function(cohort_nm) {
    if (cohort_nm == "Overall") {
      vars <- unique(unlist(default_model_vars_by_cohort[c("FE", "LTC", "YED")], use.names = FALSE))
    } else {
      vars <- default_model_vars_by_cohort[[cohort_nm]] %||% character()
    }
    
    vars <- setdiff(unique(vars), blacklist)
    vars <- intersect(vars, names(cohorts[[cohort_nm]]))
    vars
  })
  names(model_vars_by_cohort) <- cohort_names
  
  # If the caller supplies variables explicitly, restrict each cohort to that
  # requested subset after applying the blacklist.
  if (!(is.null(continuous_vars) && is.null(categorical_vars))) {
    if (is.null(continuous_vars)) continuous_vars <- character()
    if (is.null(categorical_vars)) categorical_vars <- character()
    user_vars <- setdiff(unique(c(continuous_vars, categorical_vars)), blacklist)
    model_vars_by_cohort <- lapply(model_vars_by_cohort, function(vars) intersect(vars, user_vars))
  }
  
  # Build a per-cohort lookup table saying whether each analysis variable is
  # continuous or categorical.
  type_lookup_by_cohort <- lapply(cohort_names, function(cohort_nm) {
    infer_types_for_vars(cohorts[[cohort_nm]], model_vars_by_cohort[[cohort_nm]])
  })
  names(type_lookup_by_cohort) <- cohort_names
  
  # Store the final global variable lists in the returned config for traceability.
  continuous_vars <- unique(unlist(lapply(type_lookup_by_cohort, function(x) x$variable[x$variable_type == "continuous"]), use.names = FALSE))
  categorical_vars <- unique(unlist(lapply(type_lookup_by_cohort, function(x) x$variable[x$variable_type == "categorical"]), use.names = FALSE))
  
  # ---------------------------------------------------------------------------
  # Main drift analysis
  # Iterate over cohorts and variables, compare month1 with month2, and collect
  # one row per cohort-variable pair in `drift_rows`. Category-level proportions
  # are kept separately in `props_rows`.
  # ---------------------------------------------------------------------------
  drift_rows <- list()
  props_rows <- list()
  idx <- 1L
  
  # Analyse each eligible cohort independently because each cohort uses its own
  # variable set and weight column.
  for (cohort_nm in cohort_names) {
    dfc <- cohorts[[cohort_nm]]
    wcol <- cohort_weight_map[[cohort_nm]]
    
    if (!wcol %in% names(dfc)) {
      warning(sprintf("Weight column '%s' not found for cohort '%s'; using equal weights.", wcol, cohort_nm))
    }
    
    # Add a temporary standardised weight column so downstream code does not
    # need to know the cohort-specific weight-column name.
    dfc <- dfc %>% mutate(`..w..` = get_w_from_col(., wcol))
    
    # Extract the two monthly snapshots being compared.
    d1 <- dfc %>% filter(date == month1)
    d2 <- dfc %>% filter(date == month2)
    
    if (!nrow(d1) || !nrow(d2)) {
      warning(sprintf("Skipping cohort '%s': one or both comparison months are absent.", cohort_nm))
      next
    }
    
    analysis_vars <- model_vars_by_cohort[[cohort_nm]]
    type_lookup <- type_lookup_by_cohort[[cohort_nm]]
    
    if (!length(analysis_vars)) {
      warning(sprintf("Skipping cohort '%s': no model variables available for drift analysis.", cohort_nm))
      next
    }
    
    # Compare every selected variable for this cohort.
    for (v in analysis_vars) {
      if (v %notin% names(d1) || v %notin% names(d2)) next
      
      # Missingness is recorded for both continuous and categorical variables,
      # but the drift distance itself is calculated on non-missing weighted rows.
      vtype <- type_lookup$variable_type[type_lookup$variable == v][1]
      miss1 <- mean(is.na(d1[[v]]))
      miss2 <- mean(is.na(d2[[v]]))
      
      # Continuous variables: weighted KS distance plus weighted location/spread
      # summaries for each month.
      if (identical(vtype, "continuous")) {
        dist_val <- weighted_ks_distance(d1[[v]], d1$`..w..`, d2[[v]], d2$`..w..`)
        s1 <- summarise_continuous(d1[[v]], d1$`..w..`)
        s2 <- summarise_continuous(d2[[v]], d2$`..w..`)
        
        drift_rows[[idx]] <- tibble(
          cohort = cohort_nm,
          weight_col = wcol,
          variable = v,
          variable_type = "continuous",
          metric = "weighted_ks",
          distance = dist_val,
          month1 = month1_lab,
          month2 = month2_lab,
          n_month1 = nrow(d1),
          n_month2 = nrow(d2),
          missing_prop_month1 = miss1,
          missing_prop_month2 = miss2,
          mean_w_month1 = s1$mean_w,
          mean_w_month2 = s2$mean_w,
          q25_w_month1 = s1$q25_w,
          q25_w_month2 = s2$q25_w,
          median_w_month1 = s1$median_w,
          median_w_month2 = s2$median_w,
          q75_w_month1 = s1$q75_w,
          q75_w_month2 = s2$q75_w,
          min_month1 = s1$min,
          min_month2 = s2$min,
          max_month1 = s1$max,
          max_month2 = s2$max,
          top_level_month1 = NA_character_,
          top_level_month2 = NA_character_,
          top_prop_w_month1 = NA_real_,
          top_prop_w_month2 = NA_real_
        )
      } else {
        # Categorical variables: weighted total variation distance plus top-level
        # summaries and full weighted level proportions for each month.
        dist_val <- weighted_tv_distance(d1[[v]], d1$`..w..`, d2[[v]], d2$`..w..`)
        s1 <- summarise_categorical(d1[[v]], d1$`..w..`)
        s2 <- summarise_categorical(d2[[v]], d2$`..w..`)
        
        drift_rows[[idx]] <- tibble(
          cohort = cohort_nm,
          weight_col = wcol,
          variable = v,
          variable_type = "categorical",
          metric = "weighted_tv",
          distance = dist_val,
          month1 = month1_lab,
          month2 = month2_lab,
          n_month1 = nrow(d1),
          n_month2 = nrow(d2),
          missing_prop_month1 = miss1,
          missing_prop_month2 = miss2,
          mean_w_month1 = NA_real_,
          mean_w_month2 = NA_real_,
          q25_w_month1 = NA_real_,
          q25_w_month2 = NA_real_,
          median_w_month1 = NA_real_,
          median_w_month2 = NA_real_,
          q75_w_month1 = NA_real_,
          q75_w_month2 = NA_real_,
          min_month1 = NA_real_,
          min_month2 = NA_real_,
          max_month1 = NA_real_,
          max_month2 = NA_real_,
          top_level_month1 = s1$summary$top_level,
          top_level_month2 = s2$summary$top_level,
          top_prop_w_month1 = s1$summary$top_prop_w,
          top_prop_w_month2 = s2$summary$top_prop_w
        )
        
        props_rows[[length(props_rows) + 1L]] <- bind_rows(
          s1$proportions %>% mutate(cohort = cohort_nm, variable = v, month = month1_lab),
          s2$proportions %>% mutate(cohort = cohort_nm, variable = v, month = month2_lab)
        )
      }
      
      idx <- idx + 1L
    }
  }
  
  # Combine collected rows into final tables, ordered so the largest distances
  # are easiest to inspect within each cohort/type.
  drift_tbl <- bind_rows(drift_rows) %>%
    arrange(cohort, variable_type, desc(distance), variable)
  
  cat_props_tbl <- bind_rows(props_rows) %>%
    select(cohort, variable, month, level, prop_w) %>%
    arrange(cohort, variable, month, desc(prop_w), level)
  
  # ---------------------------------------------------------------------------
  # Save tables ONLY inside out_dir
  # Keeping all outputs under `out_dir` avoids accidentally writing artefacts into
  # the project root or current working directory.
  # ---------------------------------------------------------------------------
  summary_csv <- file.path(out_dir, summary_csv_name)
  cat_props_csv <- file.path(out_dir, categorical_props_csv_name)
  
  if (isTRUE(write_outputs)) {
    # Use base write.csv for portability and UTF-8 encoding for safer handling of
    # any non-ASCII category labels.
    suppressWarnings(write.csv(drift_tbl, summary_csv, row.names = FALSE, fileEncoding = "UTF-8"))
    suppressWarnings(write.csv(cat_props_tbl, cat_props_csv, row.names = FALSE, fileEncoding = "UTF-8"))
  }
  
  # ---------------------------------------------------------------------------
  # Plots for RAW variables
  # Plot only the largest drift distances to keep output manageable. Continuous
  # variables get histogram plots; categorical variables get a weighted proportion
  # comparison plot.
  # ---------------------------------------------------------------------------
  plot_paths <- character()
  
  if (isTRUE(make_plots) && nrow(drift_tbl)) {
    # Select the top N drifting variables within each cohort and variable type.
    top_tbl <- drift_tbl %>%
      group_by(cohort, variable_type) %>%
      slice_max(order_by = distance, n = plot_top_n, with_ties = FALSE) %>%
      ungroup()
    
    # Regenerate the relevant monthly cohort subsets and save plot files.
    for (i in seq_len(nrow(top_tbl))) {
      rr <- top_tbl[i, ]
      cohort_nm <- rr$cohort
      v <- rr$variable
      dfc <- cohorts[[cohort_nm]] %>% mutate(`..w..` = get_w_from_col(., cohort_weight_map[[cohort_nm]]))
      # Extract the two monthly snapshots being compared.
      d1 <- dfc %>% filter(date == month1)
      d2 <- dfc %>% filter(date == month2)
      
      if (rr$variable_type == "continuous") {
        out_path_hist <- file.path(
          out_dir,
          "plots",
          "continuous",
          "histogram",
          paste0(cohort_nm, "__", rr$metric, "__histogram__", safe_filename(v), ".png")
        )
        plot_continuous_histogram(d1, d2, v, month1_lab, month2_lab, out_path_hist)
        plot_paths <- c(plot_paths, out_path_hist)
      } else {
        out_path <- file.path(
          out_dir,
          "plots",
          "categorical",
          paste0(cohort_nm, "__", rr$metric, "__", safe_filename(v), ".png")
        )
        plot_categorical_compare(d1, d2, v, month1_lab, month2_lab, rr$distance, out_path)
        plot_paths <- c(plot_paths, out_path)
      }
    }
  }
  
  if (isTRUE(write_outputs)) {
    # Print concise output locations for batch logs / interactive use.
    message(sprintf("[weightedDataDrift] Saved drift summary to: %s", summary_csv))
    message(sprintf("[weightedDataDrift] Saved categorical proportions to: %s", cat_props_csv))
    if (isTRUE(make_plots)) {
      message(sprintf("[weightedDataDrift] Saved plots under: %s", file.path(out_dir, "plots")))
    }
  }
  
  # Return all key artefacts invisibly so the function is convenient in scripts
  # but can still be assigned to an object for inspection/testing.
  invisible(list(
    drift_summary = drift_tbl,
    categorical_proportions = cat_props_tbl,
    plot_paths = plot_paths,
    cohorts = cohorts,
    output_files = list(
      drift_summary_csv = if (isTRUE(write_outputs)) summary_csv else NA_character_,
      categorical_props_csv = if (isTRUE(write_outputs)) cat_props_csv else NA_character_
    ),
    config = list(
      month1 = month1_lab,
      month2 = month2_lab,
      include_overall = include_overall,
      continuous_vars = continuous_vars,
      categorical_vars = categorical_vars,
      cohort_weight_map = cohort_weight_map
    )
  ))
}