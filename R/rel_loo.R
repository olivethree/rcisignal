#' Leave-one-out influence screening
#'
#' @description
#' Screens for influential producers by asking: "if this producer
#' is removed from the sample, how much does the group classification
#' image change?" Producers whose removal moves the group pattern
#' disproportionately are flagged for inspection. This is an
#' **influence / outlier screening** tool, not a reliability
#' statistic.
#'
#' @details
#' For each producer `i`, compute the Pearson correlation between
#' the full-sample group CI and the group CI recomputed without
#' that producer:
#' ```
#' full        <- rowMeans(signal_matrix)
#' r_loo[i]    <- cor(full, rowMeans(signal_matrix[, -i]))
#' ```
#'
#' Because the full-sample mean and the leave-one-out mean share
#' `(N - 1) / N` of their data, `r_loo` values are near 1 by
#' construction even on noisy data, typically `[0.95, 0.999]` at
#' `N = 30`. The absolute level of `r_loo` is not informative.
#' What is informative is the relative ordering: producers whose
#' `r_loo` sits clearly below the pack are candidates for
#' inspection. The function therefore returns a z-scored version of
#' `r_loo` in `$z_scores` using the same centre / spread estimators
#' as the flagging rule. `$z_scores` is the recommended quantity to
#' plot or report.
#'
#' Two flagging rules:
#' * `"mad"` (default): flag producers with
#'   `r_loo < median(r) - flag_threshold * mad(r)`. Robust to the
#'   few atypical producers RC datasets often contain.
#' * `"sd"`: flag producers with
#'   `r_loo < mean(r) - flag_threshold * sd(r)`. Sensitive to the
#'   very outliers it is trying to detect; emits a one-time
#'   per-session deprecation message and is scheduled for removal
#'   in v0.2.0.
#'
#' Default `flag_threshold = 2.5` is calibrated so a 30-producer
#' dataset flags roughly 0.3 producers by chance under `"sd"`,
#' rather than the ~1.5 a 2-SD rule would produce. Under `"mad"`
#' this is roughly comparable thanks to MAD's 1.4826 consistency
#' factor.
#'
#' @section What this function is, and is not:
#' `rel_loo()` is an influence-screening diagnostic. It answers
#' "which producers disproportionately shape the group CI?", not
#' "how reliable is the group CI?". For reliability use
#' [rel_split_half()] or [rel_icc()]. A flag does not mean the
#' producer is "bad"; it means the producer's individual CI sits
#' far enough from the group pattern that the data deserve a second
#' look.
#'
#' @section Reading the result:
#' * `$z_scores`, named numeric vector, per-producer standardised
#'   influence. The recommended quantity to plot or threshold.
#' * `$correlations`, named numeric vector, raw per-producer
#'   `r_loo` values.
#' * `$mean_r`, `$sd_r`, `$median_r`, `$mad_r`, centre / spread
#'   under each rule.
#' * `$threshold`, raw cutoff value on `r_loo`.
#' * `$flagged`, character vector of producer ids below threshold.
#' * `$summary_df`, one row per producer with `correlation`,
#'   `z_score`, and `flag`, sorted by `z_score`.
#' * `$flag_method`, `$flag_threshold`.
#'
#' @section Common mistakes:
#' * Reading `r_loo` as a reliability. An `r_loo` of .98 does not
#'   mean the CI is 98% reliable; it means a single producer's
#'   removal changed the group mean by 2%.
#' * Treating `$flagged` as "drop these producers". Investigate
#'   first; cross-check with the rcisignal input-side
#'   `run_diagnostics()` to rule out response-coding errors.
#' * Lowering `flag_threshold` below 2 to flag more producers; that
#'   trades real signal for noise. Use `flag_method = "mad"`
#'   instead if the SD rule is dominated by outliers.
#'
#' @section Reliability metrics expect raw masks:
#' Operates on the raw mask; results may be distorted if
#' `signal_matrix` was extracted from rendered (scaled) PNGs.
#'
#' @param signal_matrix Pixels x participants, base-subtracted.
#' @param flag_threshold Numeric multiplier on `sd` (or `mad`)
#'   below the centre. Default 2.5.
#' @param flag_method One of `"mad"` (default) or `"sd"`. SD/mean
#'   is retained for backwards compatibility and emits a one-time
#'   per-session deprecation message.
#' @param flag_threshold_sd Deprecated alias for `flag_threshold`.
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix)` (column-major) restricting the LOO
#'   correlation to a region. Build with [make_face_mask()]
#'   (parametric oval and sub-regions) or [read_face_mask()]
#'   (PNG/JPEG mask).
#' @return Object of class `rcisignal_rel_loo`.
#' @seealso [rel_loo_z()] for a tidy z-score accessor;
#'   [rel_split_half()], [rel_icc()] for reliability metrics
#'   proper; [run_reliability()].
#' @examples
#' \dontrun{
#' # In a real pipeline, signal_matrix comes from earlier steps:
#' #   signal_matrix <- ci_from_responses_briefrc(...)$signal_matrix
#' # For a self-contained demo we fabricate a small synthetic input:
#' n_pix  <- 32L * 32L
#' n_prod <- 20L
#' set.seed(1)
#' signal_matrix <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
#'
#' r <- rel_loo(signal_matrix)
#' print(r)
#' }
#' @export
rel_loo <- function(signal_matrix,
                    flag_threshold    = 2.5,
                    flag_method       = c("mad", "sd"),
                    flag_threshold_sd = NULL,
                    mask              = NULL) {
  abort_if_group_ci(signal_matrix, fn = "rel_loo")
  validate_signal_matrix(signal_matrix)
  signal_matrix <- apply_mask_to_signal(signal_matrix, mask)
  if (looks_scaled(signal_matrix)) warn_looks_scaled("signal_matrix")
  flag_method <- match.arg(flag_method)
  if (identical(flag_method, "sd")) warn_loo_sd_deprecated()
  if (!is.null(flag_threshold_sd)) {
    flag_threshold <- flag_threshold_sd
  }
  if (!is.numeric(flag_threshold) || length(flag_threshold) != 1L ||
        !is.finite(flag_threshold) || flag_threshold <= 0) {
    cli::cli_abort(
      "{.arg flag_threshold} must be a positive finite numeric scalar."
    )
  }

  n <- ncol(signal_matrix)
  if (is.null(colnames(signal_matrix))) {
    colnames(signal_matrix) <- sprintf("p%03d", seq_len(n))
  }
  full <- rowMeans(signal_matrix)
  cors <- numeric(n)
  names(cors) <- colnames(signal_matrix)
  for (i in seq_len(n)) {
    held_out <- rowMeans(signal_matrix[, -i, drop = FALSE])
    cors[i]  <- stats::cor(full, held_out)
  }

  mean_r   <- mean(cors)
  sd_r     <- stats::sd(cors)
  median_r <- stats::median(cors)
  mad_r    <- stats::mad(cors)

  z_scores <- if (flag_method == "sd") {
    if (sd_r > 0) (cors - mean_r) / sd_r else rep(0, n)
  } else {
    if (mad_r > 0) (cors - median_r) / mad_r else rep(0, n)
  }
  names(z_scores) <- names(cors)

  threshold <- if (flag_method == "sd") {
    mean_r - flag_threshold * sd_r
  } else {
    median_r - flag_threshold * mad_r
  }
  flagged <- names(cors)[cors < threshold]

  summary_df <- data.frame(
    participant_id = names(cors),
    correlation    = unname(cors),
    z_score        = unname(z_scores),
    flag           = cors < threshold,
    stringsAsFactors = FALSE
  )
  summary_df <- summary_df[order(summary_df$z_score), ,
                           drop = FALSE]
  rownames(summary_df) <- NULL

  new_rcisignal_rel_loo(
    correlations   = cors,
    z_scores       = z_scores,
    mean_r         = mean_r,
    sd_r           = sd_r,
    median_r       = median_r,
    mad_r          = mad_r,
    threshold      = threshold,
    flagged        = flagged,
    summary_df     = summary_df,
    flag_method    = flag_method,
    flag_threshold = flag_threshold
  )
}


#' Z-scored leave-one-out influence (accessor)
#'
#' @description
#' Convenience accessor that returns a data frame of producer ids
#' and their z-scored LOO influence, ordered from most-influential
#' (lowest, most negative `z_score`) to least. Accepts either a
#' signal matrix (runs [rel_loo()] under the hood) or an existing
#' `rcisignal_rel_loo` result object (cheap, no recomputation).
#'
#' @param x Either a `pixels x participants` signal matrix or an
#'   object of class `rcisignal_rel_loo` (as returned by [rel_loo()]).
#' @param ... Passed to [rel_loo()] when `x` is a signal matrix
#'   (e.g. `flag_threshold`, `flag_method`).
#' @return A data frame with columns `participant_id`,
#'   `correlation`, `z_score`, `flag`, sorted by `z_score`
#'   ascending.
#' @seealso [rel_loo()]
#' @examples
#' \dontrun{
#' n_pix  <- 32L * 32L
#' n_prod <- 20L
#' set.seed(1)
#' signal_matrix <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
#' rel_loo_z(signal_matrix)              # accept matrix directly
#' rel_loo_z(rel_loo(signal_matrix))     # or chain from a result
#' }
#' @export
rel_loo_z <- function(x, ...) {
  if (inherits(x, "rcisignal_rel_loo")) {
    return(x$summary_df)
  }
  if (is.matrix(x) && is.numeric(x)) {
    return(rel_loo(x, ...)$summary_df)
  }
  cli::cli_abort(c(
    "{.arg x} must be a numeric signal matrix or an \\
     {.cls rcisignal_rel_loo} object.",
    "i" = "Got {.cls {class(x)}}."
  ))
}
