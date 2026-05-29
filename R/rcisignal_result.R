## S3 constructors for result classes. None of the `new_*()`
## constructors are exported; users get these objects back from
## `rel_*()`, `run_*()`, `infoval()`, etc. Print / summary / plot
## methods live in R/plot_methods.R.
##
## Public constructors and predicates for the diagnostics half live
## here too: `rcisignal_diag_result()`, `rcisignal_diag_report()`,
## `is_rcisignal_diag_result()`.

#' Build an S3 result with shared metadata
#'
#' Centralised wrapper that attaches `rcisignal_version` to every
#' rcisignal result. Future fields with the same lifetime (e.g. a
#' build date or session-id) belong here.
#'
#' @keywords internal
#' @noRd
new_rcisignal_object <- function(payload, subclass) {
  structure(
    payload,
    class            = c(subclass, "rcisignal_result"),
    rcisignal_version = utils::packageVersion("rcisignal")
  )
}

# ---- Reliability (output-side) result classes -----------------------

new_rcisignal_rel_split_half <- function(r_hh, r_sb, ci_95, ci_95_sb,
                                        distribution, n_participants,
                                        n_permutations,
                                        null              = "none",
                                        null_distribution = NULL,
                                        r_hh_null_p95     = NA_real_,
                                        r_hh_excess       = NA_real_,
                                        r_sb_excess       = NA_real_) {
  new_rcisignal_object(
    list(
      r_hh              = r_hh,
      r_sb              = r_sb,
      ci_95             = ci_95,
      ci_95_sb          = ci_95_sb,
      distribution      = distribution,
      null              = null,
      null_distribution = null_distribution,
      r_hh_null_p95     = r_hh_null_p95,
      r_hh_excess       = r_hh_excess,
      r_sb_excess       = r_sb_excess,
      n_participants    = n_participants,
      n_permutations    = n_permutations
    ),
    subclass = "rcisignal_rel_split_half"
  )
}

new_rcisignal_rel_loo <- function(correlations, z_scores,
                                 mean_r, sd_r,
                                 median_r, mad_r,
                                 threshold, flagged, summary_df,
                                 flag_method, flag_threshold) {
  new_rcisignal_object(
    list(
      correlations   = correlations,
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
    ),
    subclass = "rcisignal_rel_loo"
  )
}

new_rcisignal_rel_icc <- function(icc_3_1, icc_3_k, icc_2_1, icc_2_k,
                                 ms_rows, ms_cols, ms_error,
                                 n_raters, n_targets, model, variants) {
  new_rcisignal_object(
    list(
      icc_3_1   = icc_3_1,
      icc_3_k   = icc_3_k,
      icc_2_1   = icc_2_1,
      icc_2_k   = icc_2_k,
      ms_rows   = ms_rows,
      ms_cols   = ms_cols,
      ms_error  = ms_error,
      n_raters  = n_raters,
      n_targets = n_targets,
      model     = model,
      variants  = variants
    ),
    subclass = "rcisignal_rel_icc"
  )
}

new_rcisignal_rel_cluster_test <- function(observed_t, clusters,
                                          null_distribution, img_dims,
                                          pos_labels, neg_labels,
                                          cluster_threshold, alpha,
                                          n_permutations,
                                          n_participants_a,
                                          n_participants_b,
                                          method = "threshold",
                                          paired = FALSE,
                                          tfce_map = NULL,
                                          tfce_pmap = NULL,
                                          tfce_significant_mask = NULL,
                                          tfce_H = NA_real_,
                                          tfce_E = NA_real_,
                                          tfce_n_steps = NA_integer_) {
  new_rcisignal_object(
    list(
      observed_t            = observed_t,
      method                = method,
      paired                = paired,
      clusters              = clusters,
      null_distribution     = null_distribution,
      img_dims              = img_dims,
      pos_labels            = pos_labels,
      neg_labels            = neg_labels,
      cluster_threshold     = cluster_threshold,
      alpha                 = alpha,
      n_permutations        = n_permutations,
      n_participants_a      = n_participants_a,
      n_participants_b      = n_participants_b,
      tfce_map              = tfce_map,
      tfce_pmap             = tfce_pmap,
      tfce_significant_mask = tfce_significant_mask,
      tfce_H                = tfce_H,
      tfce_E                = tfce_E,
      tfce_n_steps          = tfce_n_steps
    ),
    subclass = "rcisignal_rel_cluster_test"
  )
}

new_rcisignal_rel_dissim <- function(correlation, euclidean,
                                    euclidean_normalised,
                                    boot_cor, boot_dist,
                                    ci_cor, ci_dist,
                                    boot_se_cor, boot_se_dist,
                                    n_boot, ci_level,
                                    n_pixels,
                                    null              = "none",
                                    null_distribution = NULL,
                                    d_null_p95        = NA_real_,
                                    d_z               = NA_real_,
                                    d_ratio           = NA_real_,
                                    paired            = FALSE) {
  new_rcisignal_object(
    list(
      correlation          = correlation,
      euclidean            = euclidean,
      euclidean_normalised = euclidean_normalised,
      boot_cor             = boot_cor,
      boot_dist            = boot_dist,
      ci_cor               = ci_cor,
      ci_dist              = ci_dist,
      boot_se_cor          = boot_se_cor,
      boot_se_dist         = boot_se_dist,
      n_boot               = n_boot,
      ci_level             = ci_level,
      n_pixels             = n_pixels,
      null                 = null,
      null_distribution    = null_distribution,
      d_null_p95           = d_null_p95,
      d_z                  = d_z,
      d_ratio              = d_ratio,
      paired               = paired
    ),
    subclass = "rcisignal_rel_dissim"
  )
}

new_rcisignal_rel_infoval <- function(infoval, norms, reference,
                                     ref_median, ref_mad, trial_counts,
                                     mask, iter, n_pool, seed) {
  new_rcisignal_object(
    list(
      infoval      = infoval,
      norms        = norms,
      reference    = reference,
      ref_median   = ref_median,
      ref_mad      = ref_mad,
      trial_counts = trial_counts,
      mask         = mask,
      iter         = iter,
      n_pool       = n_pool,
      seed         = seed
    ),
    subclass = "rcisignal_rel_infoval"
  )
}

new_rcisignal_rel_report <- function(results, method, img_dims = NULL) {
  structure(
    list(
      results  = results,
      method   = method,
      img_dims = img_dims
    ),
    class            = c("rcisignal_rel_report", "rcisignal_result"),
    rcisignal_version = utils::packageVersion("rcisignal")
  )
}

# ---- Diagnostics (input-side) result classes ------------------------

#' Construct an `rcisignal_diag_result` object
#'
#' Every diagnostic check function returns an object of class
#' `"rcisignal_diag_result"`. This constructor is the single source of
#' truth for the result shape.
#'
#' @param status One of `"pass"`, `"warn"`, `"fail"`, or `"skip"`.
#' @param label Short (one-line) description of the check.
#' @param detail Human-readable explanation. Multiple lines allowed.
#' @param data Optional list of supporting data (data frames, flagged
#'   participant ids, summary statistics). Defaults to empty.
#' @return An object of class `"rcisignal_diag_result"`: a list with
#'   elements `status`, `label`, `detail`, and `data`.
#' @examples
#' rcisignal_diag_result("pass", "Response coding",
#'                      "All responses coded {-1, 1}.")
#' @export
rcisignal_diag_result <- function(status, label, detail, data = list()) {
  valid <- c("pass", "warn", "fail", "skip")
  if (!is.character(status) || length(status) != 1L ||
        !status %in% valid) {
    cli::cli_abort("{.arg status} must be one of {.val {valid}}.")
  }
  if (!is.character(label) || length(label) != 1L) {
    cli::cli_abort("{.arg label} must be a single string.")
  }
  if (!is.character(detail)) {
    cli::cli_abort("{.arg detail} must be a character vector.")
  }
  if (!is.list(data)) {
    cli::cli_abort("{.arg data} must be a list.")
  }
  structure(
    list(status = status, label = label, detail = detail, data = data),
    class = c("rcisignal_diag_result", "rcisignal_result")
  )
}

#' Test whether an object is an `rcisignal_diag_result`
#'
#' @param x Any R object.
#' @return `TRUE` if `x` inherits from `"rcisignal_diag_result"`,
#'   else `FALSE`.
#' @examples
#' r <- rcisignal_diag_result("pass", "Demo check", "Looks fine.")
#' is_rcisignal_diag_result(r)   # TRUE
#' is_rcisignal_diag_result(42)  # FALSE
#' @export
is_rcisignal_diag_result <- function(x) {
  inherits(x, "rcisignal_diag_result")
}

#' @export
print.rcisignal_diag_result <- function(x, ...) {
  cat(format_status_tag(x$status), " ", x$label, "\n", sep = "")
  for (line in x$detail) {
    cat("  ", line, "\n", sep = "")
  }
  invisible(x)
}

#' @export
format.rcisignal_diag_result <- function(x, ...) {
  c(
    paste0(format_status_tag(x$status), " ", x$label),
    paste0("  ", x$detail)
  )
}

format_status_tag <- function(status) {
  switch(
    status,
    pass = cli::col_green("[PASS]"),
    warn = cli::col_yellow("[WARN]"),
    fail = cli::col_red("[FAIL]"),
    skip = cli::col_silver("[SKIP]"),
    paste0("[", toupper(status), "]")
  )
}

#' Construct an `rcisignal_diag_report` object
#'
#' Collects the outputs of multiple diagnostic checks into a single
#' printable summary, with the method that was run and the names of
#' checks that were skipped.
#'
#' @param results A named list of [rcisignal_diag_result()] objects.
#' @param skipped_checks Character vector of check names that were
#'   not executed. Defaults to empty.
#' @param method `"2ifc"` or `"briefrc"`.
#' @return An object of class `"rcisignal_diag_report"`.
#' @examples
#' r <- rcisignal_diag_report(
#'   results = list(
#'     a = rcisignal_diag_result("pass", "Check A", "Looks fine.")
#'   ),
#'   method = "2ifc"
#' )
#' print(r)
#' @export
rcisignal_diag_report <- function(results,
                                 skipped_checks = character(),
                                 method = c("2ifc", "briefrc")) {
  method <- match.arg(method)
  if (!is.list(results) ||
        !all(vapply(results, is_rcisignal_diag_result, logical(1L)))) {
    cli::cli_abort(
      "{.arg results} must be a list of \\
       {.cls rcisignal_diag_result} objects."
    )
  }
  structure(
    list(
      results        = results,
      skipped_checks = as.character(skipped_checks),
      method         = method
    ),
    class            = c("rcisignal_diag_report", "rcisignal_result"),
    rcisignal_version = utils::packageVersion("rcisignal")
  )
}

#' @export
print.rcisignal_diag_report <- function(x, ...) {
  cat("== Data-quality report (", x$method, ") ==\n", sep = "")
  cat("\n")
  status_counts <- table(factor(
    vapply(x$results, function(r) r$status, character(1L)),
    levels = c("pass", "warn", "fail", "skip")
  ))
  for (nm in names(x$results)) {
    print(x$results[[nm]])
  }
  cat("\nSummary: ")
  cat(paste(
    sprintf("%s=%d", names(status_counts), as.integer(status_counts)),
    collapse = ", "
  ))
  cat("\n")
  if (length(x$skipped_checks) > 0L) {
    cat("\nSkipped checks:\n")
    for (nm in x$skipped_checks) {
      cat("  - ", nm, "\n", sep = "")
    }
  }
  invisible(x)
}

#' @export
summary.rcisignal_diag_report <- function(object, ...) {
  statuses <- vapply(object$results, function(r) r$status,
                     character(1L))
  data.frame(
    check  = names(object$results),
    status = statuses,
    label  = vapply(object$results, function(r) r$label,
                    character(1L)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# Diagnostic-side per-producer infoVal result. Wrapped by
# `infoval_report()` in an `rcisignal_diag_result` so it slots
# into `rcisignal_diag_report`.
new_rcisignal_diag_infoval <- function(infoval, norms, reference,
                                      ref_median, ref_mad,
                                      trial_counts, mask, iter,
                                      n_pool, seed) {
  structure(
    list(
      infoval      = infoval,
      norms        = norms,
      reference    = reference,
      ref_median   = ref_median,
      ref_mad      = ref_mad,
      trial_counts = trial_counts,
      mask         = mask,
      iter         = iter,
      n_pool       = n_pool,
      seed         = seed
    ),
    class            = c("rcisignal_diag_infoval", "rcisignal_result"),
    rcisignal_version = utils::packageVersion("rcisignal")
  )
}

#' @export
print.rcisignal_diag_infoval <- function(x, ...) {
  cat("<rcisignal_diag_infoval>\n")
  cat(sprintf("  producers : %d\n", length(x$infoval)))
  cat(sprintf("  pool size : %d\n", x$n_pool))
  cat(sprintf("  iter      : %d\n", x$iter))
  cat(sprintf(
    "  mask      : %s\n",
    if (is.null(x$mask)) "none" else
      sprintf("%d / %d pixels", sum(x$mask), length(x$mask))
  ))
  q <- stats::quantile(x$infoval, c(0, 0.25, 0.5, 0.75, 1),
                       na.rm = TRUE)
  cat(sprintf(
    "  z summary : min %.2f | Q1 %.2f | median %.2f | Q3 %.2f | max %.2f\n",
    q[1L], q[2L], q[3L], q[4L], q[5L]
  ))
  invisible(x)
}

# ---- Simulation result class ----------------------------------------

#' @keywords internal
#' @noRd
new_rcisignal_sim <- function(data, noise_matrix, base_face,
                              params, p, signal, meta,
                              rdata_path      = NULL,
                              stimuli         = NULL,
                              base_image_path = NULL) {
  structure(
    list(
      data            = data,
      noise_matrix    = noise_matrix,
      base_face       = base_face,
      params          = params,
      p               = p,
      signal          = signal,
      rdata_path      = rdata_path,
      base_image_path = base_image_path,
      stimuli         = stimuli,
      meta            = meta
    ),
    class             = c("rcisignal_sim", "rcisignal_result"),
    rcisignal_version = utils::packageVersion("rcisignal")
  )
}

#' @export
print.rcisignal_sim <- function(x, ...) {
  m <- x$meta
  cat(sprintf("<rcisignal_sim: %s>\n", m$method))
  cat(sprintf("  participants : %d  (%d per condition x %d)\n",
              m$n_per_condition * length(m$conditions),
              m$n_per_condition, length(m$conditions)))
  cat(sprintf("  conditions   : %s\n",
              paste(m$conditions, collapse = ", ")))
  cat(sprintf("  trials/p     : %d\n", m$n_trials))
  if (identical(m$method, "briefrc")) {
    cat(sprintf("  images/trial : %d  (= %d pairs)\n",
                m$images_per_trial, m$images_per_trial %/% 2L))
    cat(sprintf("  noise pool   : %d\n", m$noise_pool_size))
  } else {
    cat(sprintf("  noise pool   : %d  (one pattern per trial)\n",
                m$n_trials))
  }
  cat(sprintf("  image size   : %d x %d px\n",
              m$img_size, m$img_size))
  cat(sprintf("  signal       : %s @ %s\n",
              format(m$signal_strength), m$signal_region))
  cat(sprintf("  seed         : %d\n", m$seed))
  cat(sprintf("  elapsed      : %.1f s\n", m$elapsed_secs))
  cat(sprintf("  data         : %d rows x %d cols\n",
              nrow(x$data), ncol(x$data)))
  if (!is.null(x$rdata_path)) {
    cat(sprintf("  rdata_path   : %s\n", x$rdata_path))
  }
  if (!is.null(x$base_image_path)) {
    cat(sprintf("  base_image   : %s\n", x$base_image_path))
  }
  if (!is.null(x$stimuli)) {
    cat("  stimuli      : in-memory (portable across R sessions)\n")
  }
  invisible(x)
}

#' @export
summary.rcisignal_sim <- function(object, ...) {
  print(object, ...)
  invisible(object$data)
}
