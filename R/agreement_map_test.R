#' Per-pixel inferential agreement map with FWER control
#'
#' @description
#' Within a single condition, tests at each pixel whether the
#' producer-level signal differs from zero (one-sample t). The
#' permutation null is built by random sign-flip per producer
#' (exact under the assumption that, under H0, the producer's
#' signal contribution is symmetric around zero). Family-wise
#' error is controlled across pixels by the maximum |t| statistic.
#'
#' Use this when you want a per-pixel inferential overlay on a
#' descriptive agreement-map plot, typically paired with
#' [plot_ci_overlay()] so the significance contours are rendered
#' on top of the observed group CI.
#'
#' @section Reading the result and the paired plot:
#' * `$observed_t`: per-pixel one-sample t (sign indicates the
#'   direction of producer agreement; magnitude indicates strength).
#' * `$pmap`: per-pixel p under the max-|t| null. Already
#'   FWE-corrected at the alpha you chose; no further adjustment
#'   needed.
#' * `$significant_mask`: logical vector of pixels with
#'   `pmap < alpha`. This is the field that
#'   `plot_ci_overlay(..., test = result)` traces with black
#'   contours on top of the group CI overlay.
#' * Pass the result to `plot_ci_overlay(signal_matrix, base_image,
#'   test = result)` to get the canonical figure: blue = positive
#'   producer-mean signal, red = negative, opacity = magnitude,
#'   black contours = FWE-significant pixels at this test's alpha.
#'   To plot only the inferential mask (no continuous overlay),
#'   pass `signal_matrix = result$observed_t * result$significant_mask`
#'   or use `mask = result$significant_mask` to clip the overlay
#'   to the significant region.
#'
#' @param signal_matrix Pixels x participants, base-subtracted.
#' @param n_permutations Integer. Number of sign-flip iterations.
#'   Default 5000.
#' @param alpha Numeric. Significance level. Default 0.05.
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix)` (column-major). When supplied, the test
#'   is computed on the masked pixel subset. Pixels outside the
#'   mask are returned as `NA_real_` per-pixel and `FALSE` in the
#'   significant mask. Build with [make_face_mask()] (parametric
#'   oval and sub-regions) or [read_face_mask()] (PNG/JPEG mask).
#' @param seed Optional integer.
#' @param progress Show a `cli` progress bar.
#' @param acknowledge_scaling Logical. Forwarded to
#'   `assert_raw_signal()`.
#' @return Object of class `rcisignal_rel_agreement_map_test` with:
#' * `$observed_t`: per-pixel one-sample t.
#' * `$pmap`: per-pixel p-value under the max-|t| null.
#' * `$significant_mask`: logical, `pmap < alpha`.
#' * `$null_distribution`: numeric vector of `max_abs_t` per
#'   permutation.
#' * `$alpha`, `$n_permutations`, `$n_participants`, `$mask`.
#' @seealso [plot_ci_overlay()], [plot_agreement_map()],
#'   [rel_cluster_test()]
#' @examples
#' \dontrun{
#' # Minimal call-signature demo with a synthetic input.
#' n_pix  <- 32L * 32L
#' n_prod <- 20L
#' set.seed(1)
#' signal_matrix <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
#' res <- agreement_map_test(signal_matrix,
#'                           n_permutations = 500L, seed = 1)
#' print(res)
#' }
#'
#' \dontrun{
#' # Same function, richer input: plant a strong signal in the eye region
#' # so the FWER-controlled test has something to detect.
#' sim <- simulate_briefrc_data(
#'   n_per_condition = 20, n_trials = 60, conditions = "target",
#'   signal_region = "eyes", signal_strength = "strong", seed = 1
#' )
#' cis <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
#' agreement_map_test(cis$signal_matrix, n_permutations = 500L, seed = 1)
#' }
#'
#' \dontrun{
#' # Canonical pairing: feed the test result to plot_ci_overlay() so the
#' # FWE-significant pixels are outlined in black on top of the group CI.
#' sim   <- simulate_briefrc_data(
#'   n_per_condition = 20, n_trials = 60, conditions = "target",
#'   signal_region = "eyes", signal_strength = "strong", seed = 1
#' )
#' cis   <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
#' agree <- agreement_map_test(cis$signal_matrix,
#'                             n_permutations = 500L, seed = 1)
#' plot_ci_overlay(cis$signal_matrix, sim$base_face,
#'                 test = agree,
#'                 main = "CI overlay + FWE-significant pixels")
#' }
#' @export
agreement_map_test <- function(signal_matrix,
                               n_permutations      = 5000L,
                               alpha               = 0.05,
                               mask                = NULL,
                               seed                = NULL,
                               progress            = TRUE,
                               acknowledge_scaling = FALSE) {
  validate_signal_matrix(signal_matrix)
  assert_raw_signal(signal_matrix, acknowledge_scaling)
  n_pix_full <- nrow(signal_matrix)
  if (!is.null(mask)) {
    if (!is.logical(mask) || length(mask) != n_pix_full) {
      cli::cli_abort(c(
        "{.arg mask} must be a logical vector of length \\
         {n_pix_full}.",
        "*" = "got length {length(mask)}"
      ))
    }
    if (sum(mask) < 4L) {
      cli::cli_abort(
        "{.arg mask} selects too few pixels ({sum(mask)})."
      )
    }
    sig_use <- signal_matrix[mask, , drop = FALSE]
  } else {
    sig_use <- signal_matrix
  }

  n_pix <- nrow(sig_use)
  n_p   <- ncol(sig_use)
  n_permutations <- as.integer(n_permutations)
  if (n_permutations < 100L) {
    cli::cli_warn(
      "{.arg n_permutations} = {n_permutations} is low; \\
       FWER control will be coarse."
    )
  }

  one_sample_t <- function(mat) {
    n <- ncol(mat)
    m <- rowMeans(mat)
    v <- rowSums((mat - m)^2) / (n - 1L)
    se <- sqrt(v / n)
    t_vec <- m / se
    t_vec[!is.finite(t_vec)] <- 0
    t_vec
  }

  observed_t <- one_sample_t(sig_use)

  pid <- progress_start(n_permutations,
                        "agreement-map permutation",
                        show = progress)
  on.exit(progress_done(pid), add = TRUE)

  null_max_abs_t <- with_seed(seed, {
    out <- numeric(n_permutations)
    for (i in seq_len(n_permutations)) {
      signs <- sample(c(-1, 1), n_p, replace = TRUE)
      flipped <- sweep(sig_use, 2L, signs, `*`)
      t_perm  <- one_sample_t(flipped)
      out[i]  <- max(abs(t_perm), na.rm = TRUE)
      progress_tick(pid)
    }
    out
  })

  obs_abs <- abs(observed_t)
  pmap_use <- (1 + vapply(obs_abs, function(t_obs) {
    sum(null_max_abs_t >= t_obs)
  }, integer(1L))) / (1 + n_permutations)
  sig_use_mask <- pmap_use < alpha

  if (!is.null(mask)) {
    full_t   <- rep(NA_real_, n_pix_full)
    full_p   <- rep(NA_real_, n_pix_full)
    full_sig <- rep(FALSE,   n_pix_full)
    full_t[mask]   <- observed_t
    full_p[mask]   <- pmap_use
    full_sig[mask] <- sig_use_mask
    observed_t <- full_t
    pmap       <- full_p
    sig_mask   <- full_sig
  } else {
    pmap     <- pmap_use
    sig_mask <- sig_use_mask
  }

  structure(
    list(
      observed_t        = observed_t,
      pmap              = pmap,
      significant_mask  = sig_mask,
      null_distribution = null_max_abs_t,
      alpha             = alpha,
      n_permutations    = n_permutations,
      n_participants    = n_p,
      mask              = mask
    ),
    class            = c("rcisignal_rel_agreement_map_test",
                         "rcisignal_result"),
    rcisignal_version = utils::packageVersion("rcisignal")
  )
}

#' @export
print.rcisignal_rel_agreement_map_test <- function(x, ...) {
  warn_known_regression(x)
  cat("<rcisignal agreement map test>\n")
  cat(sprintf("  N producers:          %d\n", x$n_participants))
  cat(sprintf("  n_permutations:       %d\n", x$n_permutations))
  cat(sprintf("  alpha (FWER):         %.3f\n", x$alpha))
  n_tested <- if (is.null(x$mask)) length(x$pmap) else sum(x$mask)
  cat(sprintf("  pixels tested:        %d\n", n_tested))
  cat(sprintf("  significant pixels:   %d / %d (%.1f%%)\n",
              sum(x$significant_mask, na.rm = TRUE),
              n_tested,
              100 * sum(x$significant_mask, na.rm = TRUE) /
                n_tested))
  cat(sprintf("  observed |t| range:   [%.2f, %.2f]\n",
              min(abs(x$observed_t), na.rm = TRUE),
              max(abs(x$observed_t), na.rm = TRUE)))
  invisible(x)
}
