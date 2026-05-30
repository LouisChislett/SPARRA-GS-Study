# =============================================================================
# metrics.R
#
# Author       : Louis Chislett
#
# Purpose :
# Provide weighted model-evaluation metric helpers for binary prediction tasks.
#
#
# Notes about weighted AUC calculations
#
#   ROC-AUC is calculated from grouped prediction scores sorted in descending
#   order. For each tied score group, half of the tied negative weight is counted,
#   matching the standard tie-handling convention used by rank-based AUC.
#
#   PR-AUC is calculated as weighted average precision: the increase in recall
#   at each score threshold is multiplied by the precision at that threshold.
# =============================================================================

suppressPackageStartupMessages({
  library(tibble)
})

# ------------------------- weighted mean helper -------------------------
# Calculate a weighted mean after dropping invalid values and non-positive
# weights. Returns NA_real_ if no valid weighted observations remain.
w_mean <- function(x, w){
  # Keep only finite, non-missing values with strictly positive weights.
  ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w) & w > 0
  
  # If every row is invalid, return a numeric NA rather than failing.
  if (!any(ok)) return(NA_real_)
  
  # Weighted mean = sum(weight * value) / sum(weight).
  sum(w[ok] * x[ok]) / sum(w[ok])
}

.col_weighted_means <- function(X, w) {
  if (!length(w) || !nrow(X)) {
    return(rep(NA_real_, ncol(X)))
  }
  
  w <- as.numeric(w)
  sw <- sum(w, na.rm = TRUE)
  if (!is.finite(sw) || is.na(sw) || sw <= 0) {
    return(rep(NA_real_, ncol(X)))
  }
  
  colSums(sweep(X, 1L, w, `*`), na.rm = TRUE) / sw
}

.col_quantiles <- function(M, probs) {
  if (!length(M)) {
    return(matrix(NA_real_, nrow = 2L, ncol = 0L))
  }
  t(apply(M, 2L, quantile, probs = probs, na.rm = TRUE, names = FALSE, type = 6))
}

# ------------------------- shared metric input cleaner -------------------------
# Standardise outcome, probability, and weight vectors for all metric functions.
# This avoids repeating the same filtering and type conversion logic in each
# metric implementation.
.metric_inputs <- function(y, p, w){
  # Apply the project-level weight sanitisation helper before metric filtering.
  # NOTE: sanitize_weights() must be defined or sourced before this script runs.
  w <- sanitize_weights(w)
  
  # Keep rows with valid predicted probabilities, valid outcomes, finite weights,
  # and strictly positive weights.
  ok <- is.finite(p) & !is.na(y) & is.finite(w) & !is.na(w) & w > 0
  
  # Return consistently typed vectors restricted to valid rows only.
  list(
    y = as.integer(y[ok]),
    p = as.numeric(p[ok]),
    w = as.numeric(w[ok])
  )
}

# ------------------------- weighted losses -------------------------
# Calculate weighted binary log loss from observed outcomes, predicted
# probabilities, and observation weights.
log_loss_w <- function(y, p, w, eps=1e-12){
  # Clean and align input vectors.
  z <- .metric_inputs(y, p, w)
  
  # Clamp probabilities away from 0 and 1 to avoid log(0), then unpack weights.
  y <- z$y; p <- pmax(pmin(z$p, 1 - eps), eps); w <- z$w
  
  # Return NA if there are no valid observations after filtering.
  if (!length(y)) return(NA_real_)
  
  # Weighted negative log-likelihood for a Bernoulli outcome.
  -sum(w * (y*log(p) + (1-y)*log(1-p))) / sum(w)
}


# ------------------------- weighted grouped score summaries -------------------------
# Build grouped summaries of tied prediction scores, sorted from highest to lowest
# score. These grouped summaries are reused by the ROC-AUC and PR-AUC functions.
.weighted_score_groups_desc <- function(y, p, w){
  # Clean and align input vectors.
  z <- .metric_inputs(y, p, w)
  y <- z$y; p <- z$p; w <- z$w
  
  # AUC and PR curves require at least one positive and one negative case.
  if (!length(y) || length(unique(y)) < 2) {
    return(list(
      ok = FALSE,
      groups = tibble(),
      w_pos_total = NA_real_,
      w_neg_total = NA_real_
    ))
  }
  
  # Sort observations by predicted score from highest risk to lowest risk.
  ord <- order(p, decreasing = TRUE)
  y <- y[ord]
  p <- p[ord]
  w <- w[ord]
  
  # Identify runs of tied prediction scores after sorting.
  rr <- rle(p)
  ends <- cumsum(rr$lengths)
  starts <- c(1L, head(ends + 1L, -1L))
  
  # For each tied score group, sum the total positive-case weight.
  grp_pos <- vapply(seq_along(starts), function(i){
    idx <- starts[i]:ends[i]
    sum(w[idx][y[idx] == 1L])
  }, numeric(1))
  
  # For each tied score group, sum the total negative-case weight.
  grp_neg <- vapply(seq_along(starts), function(i){
    idx <- starts[i]:ends[i]
    sum(w[idx][y[idx] == 0L])
  }, numeric(1))
  
  # Store one row per unique prediction score, with positive and negative weight.
  groups <- tibble(
    score = rr$values,
    w_pos = grp_pos,
    w_neg = grp_neg
  )
  
  # Calculate total positive and negative weights across all valid observations.
  w_pos_total <- sum(groups$w_pos)
  w_neg_total <- sum(groups$w_neg)
  
  # Mark the grouped object as unusable if either class has zero or invalid weight.
  if (!is.finite(w_pos_total) || !is.finite(w_neg_total) ||
      w_pos_total <= 0 || w_neg_total <= 0) {
    return(list(
      ok = FALSE,
      groups = tibble(),
      w_pos_total = w_pos_total,
      w_neg_total = w_neg_total
    ))
  }
  
  # Return the grouped score table and class-weight totals.
  list(
    ok = TRUE,
    groups = groups,
    w_pos_total = w_pos_total,
    w_neg_total = w_neg_total
  )
}

# ------------------------- manual weighted ROC-AUC -------------------------
# Calculate weighted ROC-AUC manually from grouped prediction scores.
roc_auc_w <- function(y, p, w){
  # Build the grouped score representation used for rank-based AUC.
  obj <- .weighted_score_groups_desc(y, p, w)
  if (!isTRUE(obj$ok)) return(NA_real_)
  
  # Unpack grouped score table and total class weights.
  g <- obj$groups
  w_pos_total <- obj$w_pos_total
  w_neg_total <- obj$w_neg_total
  
  # For each score group, calculate the negative weight with lower scores.
  # Because scores are sorted descending, this is the negative weight appearing
  # after the current group.
  neg_after <- rev(c(0, head(cumsum(rev(g$w_neg)), -1L)))
  
  # AUC numerator: for every positive weight in a group, count all lower-scored
  # negative weight plus half the tied negative weight.
  auc_num <- sum(g$w_pos * (neg_after + 0.5 * g$w_neg))
  
  # AUC denominator: total possible weighted positive-negative pairs.
  auc_den <- w_pos_total * w_neg_total
  
  # Guard against invalid denominators.
  if (!is.finite(auc_den) || auc_den <= 0) return(NA_real_)
  
  # Return the weighted probability that a positive has a higher score than a
  # negative, with ties counted as half.
  auc_num / auc_den
}

# ------------------------- manual weighted PR-AUC (Average Precision) -------------------------
# Build a weighted precision-recall curve with one row per unique score threshold.
pr_curve_w <- function(y, p, w){
  # Build grouped score summaries. If unusable, return an empty tibble.
  obj <- .weighted_score_groups_desc(y, p, w)
  if (!isTRUE(obj$ok)) return(tibble())
  
  # Unpack grouped scores and total positive weight.
  g <- obj$groups
  w_pos_total <- obj$w_pos_total
  
  # Moving down the ranked score list, accumulate true-positive and false-positive
  # weights at each threshold.
  tp <- cumsum(g$w_pos)
  fp <- cumsum(g$w_neg)
  
  # Return threshold-level weighted recall and precision.
  tibble(
    threshold = g$score,
    tp = tp,
    fp = fp,
    recall = tp / w_pos_total,
    precision = tp / pmax(tp + fp, 1e-12)
  )
}

# Calculate weighted PR-AUC as average precision from the weighted PR curve.
pr_auc_average_precision_w <- function(y, p, w){
  # Generate the weighted precision-recall curve.
  pr <- pr_curve_w(y, p, w)
  
  # Return NA if the PR curve could not be constructed.
  if (!nrow(pr)) return(NA_real_)
  
  # Average precision sums precision at each threshold weighted by the increase
  # in recall from the previous threshold.
  rec_prev <- c(0, head(pr$recall, -1L))
  sum((pr$recall - rec_prev) * pr$precision)
}

.metric_bundle_w <- function(y, p, w, thr = 0.5){
  # Clean and align outcome, prediction, and weight vectors using the shared
  # helper, so the summary metrics are calculated on the same valid rows.
  z <- .metric_inputs(y, p, w)
  y <- z$y; p <- z$p; w <- z$w
  
  # If no valid observations remain after filtering, return a one-row metric
  # summary with zero counts and NA-valued performance metrics.
  if (!length(y)) {
    return(tibble(
      n = 0L, n_pos = 0L, n_neg = 0L,
      w_sum = 0, w_pos = 0, w_neg = 0,
      auc = NA_real_,
      pr_auc = NA_real_,
      precision = NA_real_,
      recall = NA_real_,
      accuracy = NA_real_,
      log_loss = NA_real_
    ))
  }
  
  # Convert predicted probabilities into binary class predictions using the
  # supplied threshold. The default is 0.5.
  pred <- as.integer(p >= thr)
  
  # Weighted confusion-matrix cells:
  #   TP = predicted positive and truly positive
  #   FP = predicted positive and truly negative
  #   TN = predicted negative and truly negative
  #   FN = predicted negative and truly positive
  w_tp <- sum(w * (pred == 1L & y == 1L))
  w_fp <- sum(w * (pred == 1L & y == 0L))
  w_tn <- sum(w * (pred == 0L & y == 0L))
  w_fn <- sum(w * (pred == 0L & y == 1L))
  
  # Return a single-row tibble containing unweighted counts, weighted class
  # totals, threshold-based metrics, and probability-based metrics.
  tibble(
    n         = length(y),
    n_pos     = sum(y == 1L),
    n_neg     = sum(y == 0L),
    w_sum     = sum(w),
    w_pos     = sum(w[y == 1L]),
    w_neg     = sum(w[y == 0L]),
    
    # Ranking metrics use the continuous probabilities, not the thresholded
    # predictions above.
    auc       = roc_auc_w(y, p, w),
    pr_auc    = pr_auc_average_precision_w(y, p, w),
    
    # Threshold-based metrics are calculated from the weighted confusion matrix.
    # Each denominator is checked to avoid division by zero.
    precision = if ((w_tp + w_fp) <= 0) NA_real_ else w_tp / (w_tp + w_fp),
    recall    = if ((w_tp + w_fn) <= 0) NA_real_ else w_tp / (w_tp + w_fn),
    accuracy  = if ((w_tp + w_fp + w_tn + w_fn) <= 0) NA_real_ else (w_tp + w_tn) / (w_tp + w_fp + w_tn + w_fn),
    
    # Calibration / loss metrics are calculated from the probabilities.
    log_loss  = log_loss_w(y, p, w)
  )
}

.empty_metric_boot <- function(){
  # Return an empty bootstrap-results tibble with the same metric columns used
  # for successful bootstrap replicates. This keeps downstream return types
  # consistent when no usable bootstrap samples are available.
  tibble(
    bootstrap_rep = integer(),
    auc = numeric(),
    pr_auc = numeric(),
    precision = numeric(),
    recall = numeric(),
    accuracy = numeric(),
    log_loss = numeric()
  )
}

.add_empty_metric_ci_cols <- function(est,
                                      bootstrap_method = "shared stratified bootstrap with recalculated weights") {
  # Add confidence-interval columns to a point-estimate tibble, initialised to
  # NA. These values are later replaced when usable bootstrap replicates exist.
  est %>%
    mutate(
      auc_ci_low = NA_real_, auc_ci_high = NA_real_,
      pr_auc_ci_low = NA_real_, pr_auc_ci_high = NA_real_,
      precision_ci_low = NA_real_, precision_ci_high = NA_real_,
      recall_ci_low = NA_real_, recall_ci_high = NA_real_,
      accuracy_ci_low = NA_real_, accuracy_ci_high = NA_real_,
      log_loss_ci_low = NA_real_, log_loss_ci_high = NA_real_,
      
      # Record how the bootstrap intervals were produced and how many bootstrap
      # replicates were ultimately used
      bootstrap_method = bootstrap_method,
      bootstrap_n_used = 0L
    )
}

metric_cis_from_shared_bootstrap <- function(y,
                                             p,
                                             w_point,
                                             reps,
                                             weight = c("w_ref", "w_pop"),
                                             thr = 0.5,
                                             conf.level = 0.95,
                                             bootstrap_method = NULL,
                                             return_boot = FALSE) {
  # Select which bootstrap weight vector to use from each supplied replicate.
  # match.arg() also validates that the caller used one of the accepted names.
  weight <- match.arg(weight)
  
  # If the caller did not supply a label, construct one that records the shared
  # bootstrap design and which recalculated weights were used.
  if (is.null(bootstrap_method)) {
    bootstrap_method <- sprintf(
      "shared stratified month bootstrap with recalculated %s weights",
      ifelse(weight == "w_ref", "reference-month", "population")
    )
  }
  
  # Standardise inputs for the point estimate.
  y <- as.integer(y)
  p <- as.numeric(p)
  w_point <- sanitize_weights(w_point)
  
  # Calculate point estimates on the original data, then attach empty CI columns.
  # The CI columns stay as NA if the bootstrap reps are missing or unusable.
  est <- .metric_bundle_w(y, p, w_point, thr = thr)
  out <- .add_empty_metric_ci_cols(est, bootstrap_method = bootstrap_method)
  empty_boot <- .empty_metric_boot()
  
  # Stop early if there are no bootstrap replicates, no observations, or the
  # observed outcome does not contain both binary classes.
  if (is.null(reps) || !length(reps) || !length(y) || length(unique(stats::na.omit(y))) < 2L) {
    return(if (isTRUE(return_boot)) list(summary = out, boot = empty_boot) else out)
  }
  
  # Preallocate a list to store successful bootstrap metric rows. 'used' tracks
  # how many replicates survive validation and metric calculation.
  boots <- vector("list", length(reps))
  used <- 0L
  
  # Loop over caller-supplied shared bootstrap replicates.
  for (b in seq_along(reps)) {
    rep_b <- reps[[b]]
    
    # Each replicate must provide row indices and the requested bootstrap weight
    # vector. Skip incomplete replicates rather than failing the whole function.
    if (is.null(rep_b) || is.null(rep_b$idx) || is.null(rep_b[[weight]])) next
    
    # Extract the sampled row indices and the replicate-specific weights.
    idx <- as.integer(rep_b$idx)
    wb <- sanitize_weights(rep_b[[weight]])
    
    # Keep only valid sampled rows: index must point into y/p, bootstrap weight
    # must be finite and positive, outcome must be present, and probability finite.
    ok <- !is.na(idx) & idx >= 1L & idx <= length(y) &
      is.finite(wb) & !is.na(wb) & wb > 0 &
      !is.na(y[idx]) & is.finite(p[idx])
    
    # Skip replicates that have no usable sampled observations.
    if (!any(ok)) next
    
    # Restrict this replicate to the valid sampled rows and matching weights.
    idx <- idx[ok]
    wb <- wb[ok]
    
    # Calculate metrics for this bootstrap replicate using the sampled outcomes,
    # sampled probabilities, and replicate-specific recalculated weights.
    res_b <- .metric_bundle_w(y[idx], p[idx], wb, thr = thr)
    if (!nrow(res_b)) next
    
    # Store the replicate number alongside its metric estimates.
    used <- used + 1L
    res_b$bootstrap_rep <- b
    boots[[used]] <- res_b
  }
  
  # If every replicate was skipped, return the point estimate with empty CIs.
  if (used < 1L) {
    return(if (isTRUE(return_boot)) list(summary = out, boot = empty_boot) else out)
  }
  
  # Combine successful bootstrap rows and keep only the replicate ID plus metric
  # columns needed for interval calculation and optional return.
  boot_df <- bind_rows(boots[seq_len(used)]) %>%
    select(bootstrap_rep, auc, pr_auc, precision, recall, accuracy, log_loss)
  
  # Two-sided percentile interval tail probability. For a 95% interval, alpha is
  # 0.025 and the returned percentiles are 2.5% and 97.5%.
  alpha <- (1 - conf.level) / 2
  
  # Helper to calculate a lower/upper percentile pair for one metric, after
  # removing non-finite values. Returns NA/NA if that metric has no finite values.
  qpair <- function(x){
    x <- x[is.finite(x)]
    if (!length(x)) return(c(NA_real_, NA_real_))
    as.numeric(stats::quantile(
      x,
      probs = c(alpha, 1 - alpha),
      na.rm = TRUE,
      names = FALSE,
      type = 6
    ))
  }
  
  # Calculate bootstrap percentile intervals for each metric.
  q_auc       <- qpair(boot_df$auc)
  q_pr_auc    <- qpair(boot_df$pr_auc)
  q_precision <- qpair(boot_df$precision)
  q_recall    <- qpair(boot_df$recall)
  q_accuracy  <- qpair(boot_df$accuracy)
  q_log_loss  <- qpair(boot_df$log_loss)
  
  # Fill the summary tibble with the calculated CI bounds.
  out$auc_ci_low         <- q_auc[1]
  out$auc_ci_high        <- q_auc[2]
  out$pr_auc_ci_low      <- q_pr_auc[1]
  out$pr_auc_ci_high     <- q_pr_auc[2]
  out$precision_ci_low   <- q_precision[1]
  out$precision_ci_high  <- q_precision[2]
  out$recall_ci_low      <- q_recall[1]
  out$recall_ci_high     <- q_recall[2]
  out$accuracy_ci_low    <- q_accuracy[1]
  out$accuracy_ci_high   <- q_accuracy[2]
  out$log_loss_ci_low    <- q_log_loss[1]
  out$log_loss_ci_high   <- q_log_loss[2]
  
  # Record how many bootstrap rows contributed to the intervals.
  out$bootstrap_n_used   <- nrow(boot_df)
  
  # Optionally return the replicate-level metrics as well as the summary.
  if (isTRUE(return_boot)) {
    return(list(summary = out, boot = boot_df))
  }
  
  # Default return value is the one-row summary with point estimates and CIs.
  out
}

# Compatibility alias. This no longer generates bootstrap samples.
bootstrap_metric_cis_w <- function(y, p, w, reps, weight = c("w_ref", "w_pop"),
                                   thr = 0.5, conf.level = 0.95,
                                   bootstrap_method = NULL,
                                   return_boot = FALSE, ...) {
  # Preserve the old public function name while delegating to the newer
  # implementation that consumes precomputed shared bootstrap replicates.
  metric_cis_from_shared_bootstrap(
    y = y, p = p, w_point = w, reps = reps, weight = weight,
    thr = thr, conf.level = conf.level,
    bootstrap_method = bootstrap_method, return_boot = return_boot
  )
}

# Compatibility alias. This function consumes caller-supplied shared reps only.
metrics_from_probs_w <- function(y, p, w, reps, weight = c("w_ref", "w_pop"),
                                 thr = 0.5,
                                 bootstrap_method = NULL,
                                 return_boot = FALSE, ...) {
  # Preserve another older entry point for callers that expect
  # metrics_from_probs_w(), using the default 95% confidence level.
  metric_cis_from_shared_bootstrap(
    y = y, p = p, w_point = w, reps = reps, weight = weight,
    thr = thr, conf.level = 0.95,
    bootstrap_method = bootstrap_method, return_boot = return_boot
  )
}
