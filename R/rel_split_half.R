#' Permuted split-half reliability with Spearman-Brown correction
#'
#' @description
#' Measures how stable a condition's group-level CI is across random
#' halves of its producer set. For each of `n_permutations`
#' iterations: randomly partition the producers (columns of
#' `signal_matrix`) into two roughly-equal halves, average each half
#' to get two group vectors, and compute the Pearson correlation.
#'
#' Use this to answer "would I get the same group CI if I'd run a
#' different half of my participants?".
#'
#' For odd N, one randomly-chosen producer is dropped per
#' permutation (re-drawn each iteration, not fixed) so both halves
#' have `floor(N/2)` producers.
#'
#' Reports both the mean per-permutation `r` and the Spearman-Brown
#' corrected reliability `r_SB = (2 r) / (1 + r)`. 95% CI is the
#' 2.5 / 97.5 percentile of the permutation distribution.
#'
#' @param signal_matrix Pixels x participants, base-subtracted.
#' @param n_permutations Integer number of random splits. Default
#'   2000.
#' @param null One of `"none"` (default), `"permutation"`, or
#'   `"random_responders"`.
#'   * `"none"` skips the empirical-null computation;
#'     `$null_distribution` is `NULL`.
#'   * `"permutation"`: at each iteration, generate fresh Gaussian
#'     noise per producer (no shared spatial structure); recompute
#'     `r_hh` per iteration. The empirical distribution becomes the
#'     chance baseline.
#'   * `"random_responders"`: simulate `ncol(signal_matrix)`
#'     producers responding at chance using `genMask()` machinery
#'     on the supplied `noise_matrix`. Closer to the empirical
#'     chance baseline of an actual RC experiment.
#' @param noise_matrix Pixels x pool-size numeric matrix of noise
#'   patterns. Required when `null = "random_responders"`.
#' @param seed Optional integer; if set, results are reproducible.
#'   The caller's global RNG state is restored on exit.
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix)` (column-major) restricting computation
#'   to a region of the image. Build with [make_face_mask()]
#'   (parametric oval and sub-regions) or [read_face_mask()]
#'   (PNG/JPEG mask).
#' @param progress Show a `cli` progress bar.
#' @section Reading the result:
#' * `$r_hh`, mean per-permutation Pearson r between the two halves.
#' * `$r_sb`, Spearman-Brown projected reliability of the full
#'   sample. This is usually the headline number to report.
#' * `$ci_95`, `$ci_95_sb`, percentile 95% CIs on each.
#' * `$distribution`, the full per-permutation r vector for plotting
#'   or custom CIs.
#' * `$null` (character), the null mode used.
#' * `$null_distribution`, when `null != "none"`: numeric vector of
#'   per-iteration `r_hh` under the chosen null. `NULL` otherwise.
#' * `$r_hh_null_p95`, 95th percentile of the null distribution.
#'   `NA_real_` when `null = "none"`.
#' * `$r_hh_excess`, observed `r_hh` minus the null median.
#' * `$r_sb_excess`, the same excess applied to Spearman-Brown
#'   projected `r_sb`.
#' * `$n_participants`, `$n_permutations`, metadata.
#'
#' @section Common mistakes:
#' * Reporting `$r_hh` instead of `$r_sb` when the audience cares
#'   about the full-sample CI's reliability (almost always).
#' * Running with a small `n_permutations` (< 200) to save time and
#'   then trusting the CI bounds; CIs widen with low n_perm.
#'
#' @section Reliability metrics expect raw masks:
#' Operates on the raw mask; results may be distorted if
#' `signal_matrix` was extracted from rendered (scaled) PNGs. See
#' the rcisignal output-reliability vignette chapter on raw vs.
#' rendered.
#'
#' @return An object of class `rcisignal_rel_split_half`. Fields
#'   described in **Reading the result** above.
#' @seealso [rel_loo()], [rel_icc()], [run_reliability()]
#' @references
#' Brinkman, L., Todorov, A., & Dotsch, R. (2017). Visualising
#' mental representations: A primer on noise-based reverse
#' correlation in social psychology. *European Review of Social
#' Psychology*, 28(1), 333-361.
#' \doi{10.1080/10463283.2017.1381469}
#'
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations:
#' uses in assessing rater reliability. *Psychological Bulletin*,
#' 86(2), 420-428. \doi{10.1037/0033-2909.86.2.420}
#' @export
#' @examples
#' \dontrun{
#' # Minimal call-signature demo with a synthetic input. Split-half
#' # correlation will hover near zero because the input is pure noise.
#' n_pix  <- 32L * 32L
#' n_prod <- 20L
#' set.seed(1)
#' signal_matrix <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
#' rel_split_half(signal_matrix, n_permutations = 200L, seed = 1)
#' }
#'
#' \dontrun{
#' # Same function, richer input: plant a strong signal in the eye region.
#' # Split-half correlation should now be measurably positive.
#' sim <- simulate_briefrc_data(
#'   n_per_condition = 20, n_trials = 60, conditions = "target",
#'   signal_region = "eyes", signal_strength = "strong", seed = 1
#' )
#' cis <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
#' rel_split_half(cis$signal_matrix, n_permutations = 200L, seed = 1)
#' }
rel_split_half <- function(signal_matrix,
                           n_permutations = 2000L,
                           null           = c("none", "permutation",
                                              "random_responders"),
                           noise_matrix   = NULL,
                           mask           = NULL,
                           seed           = NULL,
                           progress       = TRUE) {
  validate_signal_matrix(signal_matrix)
  signal_matrix <- apply_mask_to_signal(signal_matrix, mask)
  if (looks_scaled(signal_matrix)) warn_looks_scaled("signal_matrix")
  null <- match.arg(null)
  if (identical(null, "random_responders") && is.null(noise_matrix)) {
    cli::cli_abort(c(
      "{.code null = \"random_responders\"} requires \\
       {.arg noise_matrix}.",
      "i" = "Pass the pixels-by-pool noise matrix used at stimulus \\
             generation (the same matrix you would pass to \\
             {.fn infoval})."
    ))
  }
  if (!is.null(noise_matrix)) {
    if (!is.matrix(noise_matrix) || !is.numeric(noise_matrix)) {
      cli::cli_abort("{.arg noise_matrix} must be a numeric matrix.")
    }
    if (nrow(noise_matrix) != nrow(signal_matrix)) {
      cli::cli_abort(c(
        "Row counts of {.arg signal_matrix} and \\
         {.arg noise_matrix} must match.",
        "*" = "signal: {nrow(signal_matrix)} pixels",
        "*" = "noise:  {nrow(noise_matrix)} pixels"
      ))
    }
  }
  n_participants <- ncol(signal_matrix)
  half_size <- n_participants %/% 2L
  odd_n <- (n_participants %% 2L) == 1L
  n_permutations <- as.integer(n_permutations)
  if (n_permutations < 10L) {
    cli::cli_warn(
      "{.arg n_permutations} = {n_permutations} is very low; \\
       CI will be unstable."
    )
  }

  total_iter <- n_permutations *
    (1L + (if (null == "none") 0L else 1L))
  pid <- progress_start(total_iter, "split-half", show = progress)
  on.exit(progress_done(pid), add = TRUE)

  observed <- with_seed(seed, {
    out <- numeric(n_permutations)
    for (i in seq_len(n_permutations)) {
      idx <- sample.int(n_participants)
      if (odd_n) {
        idx <- idx[-1L]
      }
      h1 <- idx[seq_len(half_size)]
      h2 <- idx[(half_size + 1L):length(idx)]
      v1 <- rowMeans(signal_matrix[, h1, drop = FALSE])
      v2 <- rowMeans(signal_matrix[, h2, drop = FALSE])
      out[i] <- stats::cor(v1, v2)
      progress_tick(pid)
    }
    out
  })

  null_dist <- if (null == "none") {
    NULL
  } else {
    rel_split_half_null(
      n_participants = n_participants,
      n_pixels       = nrow(signal_matrix),
      null           = null,
      noise_matrix   = noise_matrix,
      n_permutations = n_permutations,
      seed           = if (is.null(seed)) NULL else seed + 1L,
      pid            = pid
    )
  }

  r_hh <- mean(observed, na.rm = TRUE)
  r_sb <- (2 * r_hh) / (1 + r_hh)
  ci_95 <- stats::quantile(observed, c(0.025, 0.975), na.rm = TRUE,
                           names = FALSE)
  obs_sb <- (2 * observed) / (1 + observed)
  ci_95_sb <- stats::quantile(obs_sb, c(0.025, 0.975), na.rm = TRUE,
                              names = FALSE)

  if (is.null(null_dist)) {
    r_hh_null_p95 <- NA_real_
    r_hh_excess   <- NA_real_
    r_sb_excess   <- NA_real_
  } else {
    r_hh_null_p95 <- unname(stats::quantile(null_dist, 0.95,
                                            na.rm = TRUE))
    r_hh_excess   <- r_hh - stats::median(null_dist, na.rm = TRUE)
    null_sb       <- (2 * null_dist) / (1 + null_dist)
    r_sb_excess   <- r_sb - stats::median(null_sb, na.rm = TRUE)
  }

  new_rcisignal_rel_split_half(
    r_hh              = r_hh,
    r_sb              = r_sb,
    ci_95             = ci_95,
    ci_95_sb          = ci_95_sb,
    distribution      = observed,
    null              = null,
    null_distribution = null_dist,
    r_hh_null_p95     = r_hh_null_p95,
    r_hh_excess       = r_hh_excess,
    r_sb_excess       = r_sb_excess,
    n_participants    = n_participants,
    n_permutations    = n_permutations
  )
}

#' Compute only the empirical-null distribution for split-half
#'
#' Standalone access to the null-distribution simulation used by
#' `rel_split_half(null = ...)`. Useful when the user wants to
#' precompute and reuse a null across multiple observed analyses
#' (e.g., when comparing across conditions with the same producer
#' count).
#'
#' @param n_participants Integer.
#' @param n_pixels Integer.
#' @param null Either `"permutation"` or `"random_responders"`.
#' @param noise_matrix Required for `"random_responders"`.
#' @param n_permutations Integer. Default 2000.
#' @param seed Optional integer.
#' @param pid Optional internal progress-bar id; end users should
#'   leave at `NULL`.
#' @return Numeric vector of length `n_permutations`.
#' @examples
#' \dontrun{
#' null_dist <- rel_split_half_null(
#'   n_participants = 20L,
#'   n_pixels       = 32L * 32L,
#'   null           = "permutation",
#'   n_permutations = 200L,
#'   seed           = 1
#' )
#' quantile(null_dist, c(0.025, 0.975))
#' }
#' @export
rel_split_half_null <- function(n_participants,
                                n_pixels,
                                null           = c("permutation",
                                                   "random_responders"),
                                noise_matrix   = NULL,
                                n_permutations = 2000L,
                                seed           = NULL,
                                pid            = NULL) {
  null <- match.arg(null)
  if (identical(null, "random_responders") && is.null(noise_matrix)) {
    cli::cli_abort(
      "{.code null = \"random_responders\"} requires \\
       {.arg noise_matrix}."
    )
  }
  half_size <- n_participants %/% 2L
  odd_n     <- (n_participants %% 2L) == 1L
  n_permutations <- as.integer(n_permutations)

  with_seed(seed, {
    out <- numeric(n_permutations)
    for (i in seq_len(n_permutations)) {
      sim <- if (null == "permutation") {
        matrix(stats::rnorm(n_pixels * n_participants),
               n_pixels, n_participants)
      } else {
        n_pool <- ncol(noise_matrix)
        sapply(seq_len(n_participants), function(j) {
          stim <- sample.int(n_pool, n_pool, replace = FALSE)
          resp <- sample(c(-1, 1), n_pool, replace = TRUE)
          as.vector((noise_matrix[, stim, drop = FALSE] %*% resp) /
                    n_pool)
        })
      }
      idx <- sample.int(n_participants)
      if (odd_n) idx <- idx[-1L]
      h1 <- idx[seq_len(half_size)]
      h2 <- idx[(half_size + 1L):length(idx)]
      v1 <- rowMeans(sim[, h1, drop = FALSE])
      v2 <- rowMeans(sim[, h2, drop = FALSE])
      out[i] <- stats::cor(v1, v2)
      progress_tick(pid)
    }
    out
  })
}
