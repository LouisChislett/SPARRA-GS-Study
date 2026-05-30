# =============================================================================
# shap.R
#
# Author       : Louis Chislett
#
# Purpose :
# Compute exact SHAP-compatible contribution summaries for the fixed odds-ratio
# SPARRA v3 models on the log-odds scale.
#
# This script:
# (1) Builds a weighted background/reference profile for each OR model
# - Uses the model definition to identify main effects and pairwise interactions
# - Calculates weighted mean log-OR contributions in the reference data
# - Stores lookup maps for interaction terms so row-level contributions can be
#   calculated efficiently
#
# (2) Computes row-level exact contribution values
# - Decomposes each model prediction into main-effect and interaction terms
# - Centres each term against the reference/background mean contribution
# - Returns contribution columns alongside row metadata for downstream summaries
#
# (3) Summarises weighted monthly contribution means with uncertainty
# - Calculates weighted mean contribution values by month and model term
# - Produces percentile confidence intervals for monthly contribution summaries
#
# (4) Compares contribution summaries between two selected months
# - Computes post-minus-pre differences in weighted mean contribution values
# =============================================================================

# =============================================================================
# Interaction lookup helpers
# =============================================================================
# Interaction effects are stored in OR tables with one row per level combination.
# These helper functions convert those tables into named lookup vectors and then
# retrieve log-OR values for observed pairs of feature values.

.make_interaction_log_lookup <- function(tbl, need) {
  # Return an empty named vector if the interaction table is unavailable.
  if (is.null(tbl) || !nrow(tbl)) return(setNames(numeric(0), character(0)))
  
  # Build a stable composite key from the two interaction variables.
  # The carriage-return separator is unlikely to appear in ordinary factor values.
  a <- as.character(tbl[[need[[1L]]]])
  b <- as.character(tbl[[need[[2L]]]])
  key <- paste(a, b, sep = "\r")
  
  # Convert ORs to safe log-ORs so missing, invalid, or non-positive ORs are
  # treated consistently by safe_log_or().
  val <- safe_log_or(tbl$OR)
  out <- stats::setNames(val, key)
  
  # Drop rows with missing names; these cannot be looked up reliably later.
  out[!is.na(names(out))]
}

.lookup_interaction_log <- function(map, a, b) {
  # Preserve zero-length input behaviour for empty data frames.
  if (!length(a)) return(numeric(0))
  
  # Recreate the same composite key used in .make_interaction_log_lookup().
  key <- paste(as.character(a), as.character(b), sep = "\r")
  out <- unname(map[key])
  
  # Missing or invalid combinations are treated as no interaction contribution.
  out[is.na(out) | !is.finite(out)] <- 0
  as.numeric(out)
}

# =============================================================================
# Background/reference contribution calculation
# =============================================================================
# Calculates the weighted reference mean contribution for
# every main effect and interaction term. Row-level contributions are centred
# against these means so that term values are interpretable relative to the
# chosen reference/background population.

prepare_logistic_shap_background <- function(ref_df, mdl, bg_w_col = NULL) {
  # Pull the model-specific OR object, main effects, interactions, and any value
  # transformations from the shared model definition function.
  def <- model_def(mdl)
  mains <- def$mains
  inter <- def$inter
  tweaks <- if (is.null(def$tweaks)) list() else def$tweaks
  OR <- def$OR
  
  # Use supplied background weights when present; otherwise fall back to equal
  # weights. Invalid weights are zeroed, with a final fallback to equal weights.
  bg_w <- sanitize_weights(ref_df[[bg_w_col]])

  if (!any(bg_w > 0)) bg_w <- rep(1, nrow(ref_df))
  
  # Initialise Weighted mean log-OR contribution for each main effect in the background - set to 0 to start with
  main_mu <- setNames(rep(0, length(mains)), mains)
  
  for (nm in mains) {
    # Skip terms that are absent from the reference data or absent from the OR object.
    if (!nm %in% names(ref_df) || is.null(OR[[nm]])) next
    
    # Apply any model-specific transformation before looking up the OR level - OR table has character keys
    key_ref <- ref_df[[nm]]
    if (!is.null(tweaks[[nm]])) key_ref <- tweaks[[nm]](key_ref)
    key_ref <- as.character(key_ref)
    
    # Store the weighted reference mean on the log-odds scale.
    main_mu[[nm]] <- w_mean(safe_log_or(lookup_or(OR[[nm]], key_ref)), bg_w)
  }
  
  # Weighted mean log-OR contribution for each interaction, plus lookup maps used
  # later when calculating row-level interaction contributions.
  inter_mu <- setNames(rep(0, length(inter)), vapply(inter, function(it) it$name, character(1)))
  inter_maps <- vector("list", length(inter))
  names(inter_maps) <- names(inter_mu)
  
  for (it in inter) {
    need <- it$need
    if (!all(need %in% names(ref_df))) next
    
    # Extract and optionally transform both variables used by the pairwise interaction.
    a_nm <- need[[1L]]
    b_nm <- need[[2L]]
    a_ref <- ref_df[[a_nm]]
    b_ref <- ref_df[[b_nm]]
    if (!is.null(tweaks[[a_nm]])) a_ref <- tweaks[[a_nm]](a_ref)
    if (!is.null(tweaks[[b_nm]])) b_ref <- tweaks[[b_nm]](b_ref)
    
    # Cache the interaction lookup map and its weighted background mean.
    map <- .make_interaction_log_lookup(it$table, need)
    inter_maps[[it$name]] <- map
    inter_mu[[it$name]] <- w_mean(.lookup_interaction_log(map, a_ref, b_ref), bg_w)
  }
  
  # The base value is the model intercept plus the weighted background means for
  # all main and interaction terms. Row-level centred contributions add to the
  # difference from this background value on the log-odds scale.
  list(
    mains = mains,
    inter = inter,
    tweaks = tweaks,
    OR = OR,
    main_mu = main_mu,
    inter_mu = inter_mu,
    inter_maps = inter_maps,
    base_value = log(as.numeric(OR$const)) + sum(main_mu, na.rm = TRUE) + sum(inter_mu, na.rm = TRUE)
  )
}

# =============================================================================
# Row-level exact contribution calculation
# =============================================================================
# For each row, this function calculates the log-OR contribution of every model
# term and subtracts the corresponding background mean. The output keeps useful
# metadata columns plus one contribution column per main effect and interaction.

compute_exact_contributions <- function(df_month, mdl, w_col, ref_df = NULL, ref_w_col = w_col) {
  if (!nrow(df_month)) return(tibble())
  
  # Load model metadata and OR lookup tables.
  def <- model_def(mdl)
  mains <- def$mains
  inter <- def$inter
  tweaks <- if (is.null(def$tweaks)) list() else def$tweaks
  OR <- def$OR
  
  
  # Calculate reference month means
  bg <- prepare_logistic_shap_background(ref_df = ref_df, mdl = mdl, bg_w_col = ref_w_col)
  
  # Keep row identifiers, dates, target, raw weighting variables, and the selected
  # analysis weight column before adding contribution columns.
  meta_keep <- intersect(
    c("id", "time", "year_month", "date", "target",
      "age_weight_raw", "decile_weight_raw", "sexM", w_col),
    names(df_month)
  )
  out <- df_month %>% dplyr::select(all_of(meta_keep))
  
  # Main-effect contribution columns: observed log-OR minus background mean log-OR.
  for (nm in mains) {
    key <- df_month[[nm]]
    if (!is.null(tweaks[[nm]])) key <- tweaks[[nm]](key)
    key <- as.character(key)
    out[[nm]] <- safe_log_or(lookup_or(OR[[nm]], key)) - bg$main_mu[[nm]]
  }
  
  # Interaction contribution columns: observed pairwise interaction log-OR minus
  # the background mean interaction contribution.
  for (it in inter) {
    need <- it$need
    if (!all(need %in% names(df_month))) next
    a_nm <- need[[1L]]
    b_nm <- need[[2L]]
    a_obs <- df_month[[a_nm]]
    b_obs <- df_month[[b_nm]]
    if (!is.null(tweaks[[a_nm]])) a_obs <- tweaks[[a_nm]](a_obs)
    if (!is.null(tweaks[[b_nm]])) b_obs <- tweaks[[b_nm]](b_obs)
    
    map <- bg$inter_maps[[it$name]]
    
    # calculate row contribution and subract background mean
    out[[it$name]] <- .lookup_interaction_log(map, a_obs, b_obs) - bg$inter_mu[[it$name]]
  }
  
  # Attach metadata describing the decomposition. These attributes are useful for
  # downstream checks but do not affect the returned tibble columns.
  attr(out, "shap_base_value") <- bg$base_value
  attr(out, "shap_scale") <- "log-odds"
  attr(out, "shap_level") <- "model_term"
  out
}

# =============================================================================
# Bootstrap confidence intervals for contribution means
# =============================================================================
# This function does not resample rows itself. It consumes shared bootstrap
# replicate objects produced by the wider modelling pipeline, ensuring that SHAP
# summaries use the same bootstrap path as the model metrics.

shap_cis_from_shared_bootstrap <- function(X,
                                           w_point,
                                           reps,
                                           weight = c("w_ref", "w_pop"),
                                           conf = 0.95) {
  weight <- match.arg(weight)
  
  # Ensure contribution values are in a numeric matrix with one column per term.
  if (!is.matrix(X)) X <- as.matrix(X)
  storage.mode(X) <- "double"
  
  # Point estimate uses the observed rows and their point-estimate weights.
  est <- .col_weighted_means(X, sanitize_weights(w_point))
  ci_lo <- rep(NA_real_, ncol(X))
  ci_hi <- rep(NA_real_, ncol(X))
  
  # With no valid bootstrap replicates, return point estimates and missing CIs.
  if (is.null(reps) || !length(reps) || !nrow(X)) {
    return(list(
      estimate = est,
      ci_low = ci_lo,
      ci_high = ci_hi,
      bootstrap_n_used = 0L,
      boot = matrix(NA_real_, nrow = 0L, ncol = ncol(X))
    ))
  }
  
  # Each row of boot_mat will hold one bootstrap replicate's weighted mean for
  # every contribution column.
  boot_mat <- matrix(NA_real_, nrow = length(reps), ncol = ncol(X))
  used <- 0L
  
  for (b in seq_along(reps)) {
    rep_b <- reps[[b]]
    if (is.null(rep_b) || is.null(rep_b$idx) || is.null(rep_b[[weight]])) next
    
    idx <- as.integer(rep_b$idx)
    wb <- sanitize_weights(rep_b[[weight]])
    
    # Keep only sampled row indices that point into X and have positive finite weights.
    ok <- !is.na(idx) & idx >= 1L & idx <= nrow(X) &
      is.finite(wb) & !is.na(wb) & wb > 0
    
    if (!any(ok)) next
    
    used <- used + 1L
    boot_mat[used, ] <- .col_weighted_means(X[idx[ok], , drop = FALSE], wb[ok])
  }
  
  # If every replicate was unusable, return point estimates and missing CIs.
  if (used < 1L) {
    return(list(
      estimate = est,
      ci_low = ci_lo,
      ci_high = ci_hi,
      bootstrap_n_used = 0L,
      boot = matrix(NA_real_, nrow = 0L, ncol = ncol(X))
    ))
  }
  
  # Drop unused preallocated rows and compute percentile confidence limits.
  boot_mat <- boot_mat[seq_len(used), , drop = FALSE]
  alpha <- (1 - conf) / 2
  qs <- .col_quantiles(boot_mat, probs = c(alpha, 1 - alpha))
  
  list(
    estimate = est,
    ci_low = qs[, 1L],
    ci_high = qs[, 2L],
    bootstrap_n_used = used,
    boot = boot_mat
  )
}

# =============================================================================
# Monthly weighted SHAP/contribution summaries
# =============================================================================
# Takes row-level contribution output, splits it by month, and returns one row per
# model term per month with weighted mean contribution and bootstrap CI columns.

summarise_shap_weighted <- function(shap_df,
                                    w_col,
                                    mdl,
                                    ref_df = NULL,
                                    ref_month = "2012-01",
                                    conf = 0.95,
                                    shared_reps = NULL,
                                    weight = c("w_ref", "w_pop"),
                                    bootstrap_reweights = NULL,
                                    ...) {
  weight <- match.arg(weight)
  if (!nrow(shap_df)) return(tibble())
  
  
  # Contribution columns are all non-metadata columns. Each becomes a variable in
  # the long monthly summary table.
  value_cols <- setdiff(
    names(shap_df),
    c("id", "time", "year_month", "date", "target",
      "age_weight_raw", "decile_weight_raw", "sexM", "model", w_col)
  )
  if (!length(value_cols)) return(tibble())
  
  # Backwards-compatible name for callers that already pass a month-keyed list.
  if (is.null(shared_reps) && !is.null(bootstrap_reweights)) {
    shared_reps <- bootstrap_reweights
  }
  
  # Work month-by-month so each month uses its own shared bootstrap replicate list.
  month_splits <- split(shap_df, shap_df$year_month)
  out_parts <- vector("list", length(month_splits))
  nm_i <- 0L
  
  for (ym in names(month_splits)) {
    dm <- month_splits[[ym]]
    if (!nrow(dm)) next
    
    # Point-estimate weights and date label for this month.
    w0 <- sanitize_weights(dm[[w_col]])
    d0 <- dm$date[[1]]
    
    # Matrix of row-level contribution values for this month.
    X <- as.matrix(dm[, value_cols, drop = FALSE])
    storage.mode(X) <- "double"
    
    # Pull the matching bootstrap replicates for this month when available.
    reps_ym <- NULL
    if (!is.null(shared_reps) && !is.null(shared_reps[[ym]])) {
      reps_ym <- shared_reps[[ym]]
    }
    
    # Calculate weighted means and percentile CIs for all contribution columns.
    ci_obj <- shap_cis_from_shared_bootstrap(
      X = X,
      w_point = w0,
      reps = reps_ym,
      weight = weight,
      conf = conf
    )
    
    # Store one long-format row per contribution variable for this month.
    nm_i <- nm_i + 1L
    out_parts[[nm_i]] <- tibble(
      model = mdl,
      year_month = ym,
      date = d0,
      variable = value_cols,
      n = nrow(dm),
      w_sum = sum(w0, na.rm = TRUE),
      mean_phi = ci_obj$estimate,
      mean_phi_ci_low = ci_obj$ci_low,
      mean_phi_ci_high = ci_obj$ci_high,
      bootstrap_n_used = ci_obj$bootstrap_n_used,
      ci_conf_level = conf,
      bootstrap_method = sprintf(
        "shared stratified month bootstrap with recalculated %s weights",
        ifelse(weight == "w_ref", "reference-month", "population")
      )
    )
  }
  
  bind_rows(out_parts[seq_len(nm_i)])
}

# =============================================================================
# Pre/post weighted SHAP/contribution difference summaries
# =============================================================================
# Compares two selected months by subtracting the weighted mean contribution in
# pre_month from the weighted mean contribution in post_month. Bootstrap CIs are
# calculated from paired shared bootstrap replicates when both months are present.

compute_shap_difference_weighted <- function(shap_df,
                                             w_col,
                                             mdl,
                                             ref_df = NULL,
                                             ref_month = "2012-01",
                                             pre_month = "2020-03",
                                             post_month = "2021-03",
                                             conf = 0.95,
                                             shared_reps = NULL,
                                             weight = c("w_ref", "w_pop"),
                                             bootstrap_reweights = NULL,
                                             ...) {
  weight <- match.arg(weight)
  
  if (!nrow(shap_df)) {
    return(list(summary = tibble(), boot = tibble()))
  }
  
  
  # Contribution columns are all non-metadata columns.
  value_cols <- setdiff(
    names(shap_df),
    c("id", "time", "year_month", "date", "target",
      "age_weight_raw", "decile_weight_raw", "sexM", "model", w_col)
  )
  if (!length(value_cols)) {
    return(list(summary = tibble(), boot = tibble()))
  }
  
  # Extract the two months to compare. If either is missing, there is no valid
  # difference to calculate.
  dm_pre  <- shap_df %>% filter(year_month == pre_month)
  dm_post <- shap_df %>% filter(year_month == post_month)
  if (!nrow(dm_pre) || !nrow(dm_post)) {
    return(list(summary = tibble(), boot = tibble()))
  }
  
  # Backwards-compatible name for callers that already pass bootstrap_reweights.
  if (is.null(shared_reps) && !is.null(bootstrap_reweights)) {
    shared_reps <- bootstrap_reweights
  }
  
  # Contribution matrices for the two months.
  X_pre <- as.matrix(dm_pre[, value_cols, drop = FALSE])
  X_post <- as.matrix(dm_post[, value_cols, drop = FALSE])
  storage.mode(X_pre) <- "double"
  storage.mode(X_post) <- "double"
  
  # Point-estimate weighted means and observed post-minus-pre difference.
  w_pre <- sanitize_weights(dm_pre[[w_col]])
  w_post <- sanitize_weights(dm_post[[w_col]])
  
  mean_pre <- .col_weighted_means(X_pre, w_pre)
  mean_post <- .col_weighted_means(X_post, w_post)
  diff_obs <- mean_post - mean_pre
  
  # Retrieve the bootstrap replicate lists for the two selected months.
  reps_pre <- if (!is.null(shared_reps) && !is.null(shared_reps[[pre_month]])) shared_reps[[pre_month]] else NULL
  reps_post <- if (!is.null(shared_reps) && !is.null(shared_reps[[post_month]])) shared_reps[[post_month]] else NULL
  
  alpha <- (1 - conf) / 2
  boot_df <- tibble()
  ci_low <- rep(NA_real_, length(value_cols))
  ci_high <- rep(NA_real_, length(value_cols))
  
  # Use paired replicate positions where both months have shared bootstrap objects.
  n_boot_used <- min(length(reps_pre), length(reps_post))
  used <- 0L
  
  if (n_boot_used >= 1L) {
    boot_mat <- matrix(NA_real_, nrow = n_boot_used, ncol = length(value_cols))
    
    for (b in seq_len(n_boot_used)) {
      rep_pre <- reps_pre[[b]]
      rep_post <- reps_post[[b]]
      
      if (is.null(rep_pre$idx) || is.null(rep_pre[[weight]]) ||
          is.null(rep_post$idx) || is.null(rep_post[[weight]])) next
      
      idx_pre <- as.integer(rep_pre$idx)
      idx_post <- as.integer(rep_post$idx)
      wb_pre <- sanitize_weights(rep_pre[[weight]])
      wb_post <- sanitize_weights(rep_post[[weight]])
      
      # Keep only valid sampled rows with positive finite bootstrap weights.
      ok_pre <- !is.na(idx_pre) & idx_pre >= 1L & idx_pre <= nrow(X_pre) &
        is.finite(wb_pre) & !is.na(wb_pre) & wb_pre > 0
      ok_post <- !is.na(idx_post) & idx_post >= 1L & idx_post <= nrow(X_post) &
        is.finite(wb_post) & !is.na(wb_post) & wb_post > 0
      
      if (!any(ok_pre) || !any(ok_post)) next
      
      # Compute this replicate's post-minus-pre contribution difference.
      m_pre <- .col_weighted_means(X_pre[idx_pre[ok_pre], , drop = FALSE], wb_pre[ok_pre])
      m_post <- .col_weighted_means(X_post[idx_post[ok_post], , drop = FALSE], wb_post[ok_post])
      
      used <- used + 1L
      boot_mat[used, ] <- m_post - m_pre
    }
    
    if (used >= 1L) {
      # Drop unused rows, compute percentile intervals, and keep the replicate-
      # level differences for optional downstream inspection.
      boot_mat <- boot_mat[seq_len(used), , drop = FALSE]
      qs <- .col_quantiles(boot_mat, probs = c(alpha, 1 - alpha))
      ci_low <- qs[, 1L]
      ci_high <- qs[, 2L]
      
      boot_df <- as_tibble(boot_mat)
      names(boot_df) <- value_cols
      boot_df <- boot_df %>%
        mutate(
          model = mdl,
          pre_month = pre_month,
          post_month = post_month,
          bootstrap_rep = row_number(),
          .before = 1L
        )
    }
  }
  
  # Long-format summary table: one row per contribution variable.
  summary <- tibble(
    model = mdl,
    pre_month = pre_month,
    post_month = post_month,
    variable = value_cols,
    mean_pre = mean_pre,
    mean_post = mean_post,
    diff_post_minus_pre = diff_obs,
    diff_ci_low = ci_low,
    diff_ci_high = ci_high,
    bootstrap_n_used = used,
    ci_conf_level = conf,
    bootstrap_method = sprintf(
      "paired shared stratified month bootstrap with recalculated %s weights",
      ifelse(weight == "w_ref", "reference-month", "population")
    )
  )
  
  list(summary = summary, boot = boot_df)
}
