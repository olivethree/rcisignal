#' Between-condition dissimilarity with bootstrap confidence intervals
#'
#' @description
#' Quantifies overall dissimilarity between two conditions'
#' group-level classification images. The primary statistic is
#' Euclidean distance between the two group-mean CIs, reported
#' both raw and normalised by `sqrt(n_pixels)` for cross-resolution
#' comparability. Percentile bootstrap 95% confidence intervals are
#' computed by resampling participants with replacement within each
#' condition.
#'
#' Pair with [rel_cluster_test()] when you also want to know
#' **where** the two conditions differ.
#'
#' @details
#' For the observed statistics and each bootstrap replicate `i`:
#' ```
#' mean_a      = rowMeans(signal_matrix_a)
#' mean_b      = rowMeans(signal_matrix_b)
#' observed_dist            = sqrt(sum((mean_a - mean_b)^2))
#' observed_dist_normalised = observed_dist / sqrt(n_pixels)
#' ```
#' Percentile CI via base R `quantile()`; no `boot` dependency.
#'
#' @section Why Euclidean and not Pearson correlation as the primary:
#' Two base-subtracted CIs share systematic image-domain spatial
#' structure (face shape, oval signal support, low-frequency
#' Gaussian-noise smoothness) that pushes their correlation above
#' zero even when the underlying mental representations are
#' unrelated. Absolute correlation values therefore do not cleanly
#' mean "these conditions are similar": two unrelated traits can
#' easily produce `r` well above zero from the shared image
#' scaffolding alone, so a "high" correlation does not by itself
#' license a similarity claim. Euclidean distance does not share
#' this baseline issue, which is why it is the primary statistic
#' here.
#'
#' The Pearson correlation fields (`$correlation`, `$boot_cor`,
#' `$ci_cor`, `$boot_se_cor`) are returned alongside it as a
#' secondary summary, for users who need to report `r` for
#' comparability with prior literature or with another analysis
#' pipeline. They are not recommended as a standalone similarity
#' score. If you do report `r`, interpret it **relatively** (the
#' ordering of `r` across multiple condition pairs is more
#' defensible than any single absolute value) and **benchmark
#' against a permutation null** rather than against zero. The
#' image-domain scaffolding shifts the chance baseline upward;
#' "above zero" is not the right comparison.
#'
#' @param signal_matrix_a,signal_matrix_b Pixels x participants,
#'   base-subtracted. Row counts must match.
#' @param paired Logical. `FALSE` (default) for between-subjects:
#'   participants resampled within A and B independently. `TRUE`
#'   for within-subjects: A and B share a single resample index per
#'   replicate so the paired covariance structure is preserved.
#' @param n_boot Bootstrap replicates. Default 2000.
#' @param ci_level Confidence level. Default 0.95.
#' @param null One of `"none"` (default) or `"permutation"`.
#'   `"permutation"` builds an empirical chance baseline for the
#'   Euclidean distance.
#' @param n_permutations Integer. Number of null iterations when
#'   `null = "permutation"`. Default 2000.
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix_a)` (column-major) restricting the
#'   Euclidean / correlation computation to a region. Build with
#'   [make_face_mask()] (parametric oval and sub-regions) or
#'   [read_face_mask()] (PNG/JPEG mask).
#' @param seed Optional integer; RNG state restored on exit.
#' @param progress Show a `cli` progress bar.
#' @param acknowledge_scaling Logical. When `FALSE` (default), the
#'   shared `assert_raw_signal()` helper errors on a known-rendered
#'   matrix on either side.
#' @section Reading the plot:
#' `plot()` on the returned object renders the bootstrap
#' distributions as two side-by-side histograms (Euclidean
#' distance on the left, Pearson r on the right). The Pearson
#' panel is rendered in grey because it is a secondary summary,
#' not recommended as a standalone similarity score; see the
#' description for why and how to use it carefully if needed.
#' * The shaded vertical band marks the percentile CI at
#'   `ci_level`.
#' * The vertical line marks the *observed* statistic on the
#'   real data (not a bootstrap mean).
#' * A non-overlapping CI band away from zero on the Euclidean
#'   panel indicates the two group-mean CIs sit a non-trivial
#'   distance apart in pixel space, robust to participant-level
#'   resampling. The numbers are returned in `$ci_dist`.
#' * For visual comparison across multiple contrasts, pass each
#'   `rel_dissimilarity()` result to [plot_dissimilarity_grid()],
#'   which lays them out as labelled CI bars on a shared axis.
#' * For a *spatial* picture of where the two conditions differ,
#'   pair this with [rel_cluster_test()] (or use
#'   [run_discriminability()] to run both in one call).
#'
#' @section Reading the result:
#' * `$euclidean`, observed Euclidean distance between group means
#'   (primary statistic).
#' * `$euclidean_normalised`, `$euclidean / sqrt(n_pixels)`. Use
#'   for cross-resolution comparisons.
#' * `$boot_dist`, `$ci_dist`, `$boot_se_dist`, bootstrap
#'   distribution, percentile CI, and SE of the Euclidean distance.
#' * `$null` (character), the null mode used.
#' * `$null_distribution`, when `null != "none"`: numeric vector
#'   of per-iteration Euclidean distances under the chosen null.
#' * `$d_null_p95`, 95th percentile of the null distribution.
#' * `$d_z`, z-equivalent effect size:
#'   `(observed_d - mean(null)) / sd(null)`.
#' * `$d_ratio`, observed Euclidean over the null median.
#' * `$correlation`, `$boot_cor`, `$ci_cor`, `$boot_se_cor`:
#'   Pearson correlation of the group means and its bootstrap
#'   summaries. Secondary; not recommended as a standalone
#'   similarity score (see the "Why Euclidean" section). If
#'   reporting, prefer relative comparisons across pairs against
#'   a permutation null.
#' * `$n_boot`, `$ci_level`, `$paired`, metadata.
#'
#' @return Object of class `rcisignal_rel_dissim`.
#' @seealso [rel_cluster_test()], [run_discriminability()]
#' @references
#' Efron, B., & Tibshirani, R. J. (1994). *An introduction to the
#' bootstrap*. Chapman & Hall / CRC.
#' @examples
#' \dontrun{
#' # Minimal call-signature demo with two synthetic inputs.
#' n_pix  <- 32L * 32L
#' n_prod <- 20L
#' set.seed(1)
#' signal_matrix_a <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
#' signal_matrix_b <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
#' rel_dissimilarity(signal_matrix_a, signal_matrix_b,
#'                   n_boot = 200L, seed = 1)
#' }
#'
#' \dontrun{
#' # Same function, richer input: signal planted in different face regions
#' # (eyes vs mouth). The Euclidean distance and its bootstrap CI should
#' # be well above zero, reflecting genuine spatial divergence.
#' sim_eyes  <- simulate_briefrc_data(
#'   n_per_condition = 20, n_trials = 60, conditions = "x",
#'   signal_region = "eyes", signal_strength = "strong", seed = 1
#' )
#' sim_mouth <- simulate_briefrc_data(
#'   n_per_condition = 20, n_trials = 60, conditions = "x",
#'   signal_region = "mouth", signal_strength = "strong", seed = 2
#' )
#' sig_eyes  <- ci_from_responses_briefrc(
#'   sim_eyes$data, noise_matrix = sim_eyes$noise_matrix)$signal_matrix
#' sig_mouth <- ci_from_responses_briefrc(
#'   sim_mouth$data, noise_matrix = sim_mouth$noise_matrix)$signal_matrix
#' d <- rel_dissimilarity(sig_eyes, sig_mouth, n_boot = 500L, seed = 1)
#' # Bootstrap distribution + observed Euclidean + 95% CI band.
#' plot(d, main = "Eyes vs Mouth: bootstrap dissimilarity")
#' }
#' @export
rel_dissimilarity <- function(signal_matrix_a,
                              signal_matrix_b,
                              paired              = FALSE,
                              n_boot              = 2000L,
                              ci_level            = 0.95,
                              null                = c("none",
                                                       "permutation"),
                              n_permutations      = 2000L,
                              mask                = NULL,
                              seed                = NULL,
                              progress            = TRUE,
                              acknowledge_scaling = FALSE) {
  abort_if_group_ci(signal_matrix_a, fn = "rel_dissimilarity",
                    arg = "signal_matrix_a")
  abort_if_group_ci(signal_matrix_b, fn = "rel_dissimilarity",
                    arg = "signal_matrix_b")
  validate_two_signal_matrices(signal_matrix_a, signal_matrix_b)
  if (isTRUE(paired)) {
    validate_paired_matrices(signal_matrix_a, signal_matrix_b)
  }
  assert_raw_signal(signal_matrix_a, acknowledge_scaling,
                    name = "signal_matrix_a")
  assert_raw_signal(signal_matrix_b, acknowledge_scaling,
                    name = "signal_matrix_b")
  signal_matrix_a <- apply_mask_to_signal(signal_matrix_a, mask,
                                          name = "signal_matrix_a")
  signal_matrix_b <- apply_mask_to_signal(signal_matrix_b, mask,
                                          name = "signal_matrix_b")
  null <- match.arg(null)
  n_a      <- ncol(signal_matrix_a)
  n_b      <- ncol(signal_matrix_b)
  n_pixels <- nrow(signal_matrix_a)
  n_boot   <- as.integer(n_boot)
  n_permutations <- as.integer(n_permutations)

  mean_a <- rowMeans(signal_matrix_a)
  mean_b <- rowMeans(signal_matrix_b)

  observed_cor       <- stats::cor(mean_a, mean_b)
  observed_dist      <- sqrt(sum((mean_a - mean_b)^2))
  observed_dist_norm <- observed_dist / sqrt(n_pixels)

  boot_cor  <- numeric(n_boot)
  boot_dist <- numeric(n_boot)

  total_iter <- n_boot +
    (if (null == "none") 0L else n_permutations)
  pid <- progress_start(total_iter, "dissimilarity bootstrap",
                        show = progress)
  on.exit(progress_done(pid), add = TRUE)

  with_seed(seed, {
    for (i in seq_len(n_boot)) {
      if (isTRUE(paired)) {
        idx   <- sample.int(n_a, n_a, replace = TRUE)
        idx_a <- idx
        idx_b <- idx
      } else {
        idx_a <- sample.int(n_a, n_a, replace = TRUE)
        idx_b <- sample.int(n_b, n_b, replace = TRUE)
      }
      m_a <- rowMeans(signal_matrix_a[, idx_a, drop = FALSE])
      m_b <- rowMeans(signal_matrix_b[, idx_b, drop = FALSE])
      boot_cor[i]  <- stats::cor(m_a, m_b)
      boot_dist[i] <- sqrt(sum((m_a - m_b)^2))
      progress_tick(pid)
    }
  })

  d_null <- if (null == "none") {
    NULL
  } else {
    if (isTRUE(paired)) {
      diff_mat <- signal_matrix_a - signal_matrix_b
      with_seed(if (is.null(seed)) NULL else seed + 1L, {
        out <- numeric(n_permutations)
        for (i in seq_len(n_permutations)) {
          signs <- sample(c(-1, 1), n_a, replace = TRUE)
          d_flipped <- sweep(diff_mat, 2L, signs, `*`)
          m_d <- rowMeans(d_flipped)
          out[i] <- sqrt(sum(m_d * m_d))
          progress_tick(pid)
        }
        out
      })
    } else {
      combined <- cbind(signal_matrix_a, signal_matrix_b)
      n_total <- n_a + n_b
      with_seed(if (is.null(seed)) NULL else seed + 1L, {
        out <- numeric(n_permutations)
        for (i in seq_len(n_permutations)) {
          perm <- sample.int(n_total)
          idx_a <- perm[seq_len(n_a)]
          idx_b <- perm[(n_a + 1L):n_total]
          m_a   <- rowMeans(combined[, idx_a, drop = FALSE])
          m_b   <- rowMeans(combined[, idx_b, drop = FALSE])
          out[i] <- sqrt(sum((m_a - m_b)^2))
          progress_tick(pid)
        }
        out
      })
    }
  }

  tail <- (1 - ci_level) / 2
  ci_cor  <- stats::quantile(boot_cor,  c(tail, 1 - tail),
                             na.rm = TRUE, names = FALSE)
  ci_dist <- stats::quantile(boot_dist, c(tail, 1 - tail),
                             na.rm = TRUE, names = FALSE)

  if (is.null(d_null)) {
    d_null_p95 <- NA_real_
    d_z        <- NA_real_
    d_ratio    <- NA_real_
  } else {
    sd_null    <- stats::sd(d_null, na.rm = TRUE)
    mean_null  <- mean(d_null, na.rm = TRUE)
    d_null_p95 <- unname(stats::quantile(d_null, 0.95,
                                         na.rm = TRUE))
    d_z        <- if (sd_null > 0) {
      (observed_dist - mean_null) / sd_null
    } else NA_real_
    med_null   <- stats::median(d_null, na.rm = TRUE)
    d_ratio    <- if (med_null > 0)
      observed_dist / med_null else NA_real_
  }

  new_rcisignal_rel_dissim(
    correlation          = observed_cor,
    euclidean            = observed_dist,
    euclidean_normalised = observed_dist_norm,
    boot_cor             = boot_cor,
    boot_dist            = boot_dist,
    ci_cor               = ci_cor,
    ci_dist              = ci_dist,
    boot_se_cor          = stats::sd(boot_cor,  na.rm = TRUE),
    boot_se_dist         = stats::sd(boot_dist, na.rm = TRUE),
    n_boot               = n_boot,
    ci_level             = ci_level,
    n_pixels             = n_pixels,
    null                 = null,
    null_distribution    = d_null,
    d_null_p95           = d_null_p95,
    d_z                  = d_z,
    d_ratio              = d_ratio,
    paired               = isTRUE(paired)
  )
}
