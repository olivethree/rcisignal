#' Vectorised pixel-wise t-test (independent or paired)
#'
#' @description
#' At every pixel, tests whether condition A's mean signal differs
#' from condition B's. Two modes:
#'
#' * `paired = FALSE` (default): independent-samples Welch t per
#'   pixel. Correct when producers in A and B are different people
#'   (between-subjects design).
#' * `paired = TRUE`: paired t per pixel on the per-producer
#'   difference `A - B`. Correct when the same producers
#'   contributed to both conditions (within-subjects design).
#'
#' Returns a numeric vector of t-values, length `n_pixels`. Used as
#' an intermediate by [rel_cluster_test()]; not intended as a
#' standalone inferential test.
#'
#' @param signal_matrix_a,signal_matrix_b Pixels x participants,
#'   base-subtracted. Row counts must match. When `paired = TRUE`
#'   the column counts must also match, and column names must
#'   correspond to the same producer across matrices.
#' @param paired Logical. `FALSE` (default) uses independent
#'   Welch t; `TRUE` uses paired t.
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix_a)` (column-major). Both matrices are
#'   subsetted with the same mask before computing t; the returned
#'   vector is then of length `sum(mask)`. Build with
#'   [make_face_mask()] (parametric oval and sub-regions) or
#'   [read_face_mask()] (PNG/JPEG mask).
#' @param acknowledge_scaling Logical. When `FALSE` (default), the
#'   shared `assert_raw_signal()` helper errors on a known-rendered
#'   matrix. Cascades to internal `pixel_t_test()` calls from
#'   [rel_cluster_test()].
#' @return Numeric vector of length `nrow(signal_matrix_a)` (or
#'   `sum(mask)` if `mask` is supplied). Pixels with zero variance
#'   get `0` rather than `NaN`.
#' @section Reliability metrics expect raw masks:
#' Welch t and paired t are variance-based and sensitive to
#' scaling. Inputs with `attr(., "source") == "rendered"` (set
#' automatically by Mode 1 readers like [extract_signal()]) error
#' unless `acknowledge_scaling = TRUE`.
#' @seealso [rel_cluster_test()]
#' @examples
#' \dontrun{
#' # Minimal call-signature demo with two synthetic inputs.
#' n_pix  <- 32L * 32L
#' n_prod <- 20L
#' set.seed(1)
#' signal_matrix_a <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
#' signal_matrix_b <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
#' summary(pixel_t_test(signal_matrix_a, signal_matrix_b))
#' }
#'
#' \dontrun{
#' # Same function, richer input: simulate two conditions with signal
#' # planted in different face regions (eyes vs mouth). The pixel-wise
#' # t-map should be large-positive around the eyes and large-negative
#' # around the mouth — pixels where the conditions disagree.
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
#' summary(pixel_t_test(sig_eyes, sig_mouth))
#' }
#' @export
pixel_t_test <- function(signal_matrix_a, signal_matrix_b,
                         paired              = FALSE,
                         mask                = NULL,
                         acknowledge_scaling = FALSE) {
  abort_if_group_ci(signal_matrix_a, fn = "pixel_t_test",
                    arg = "signal_matrix_a")
  abort_if_group_ci(signal_matrix_b, fn = "pixel_t_test",
                    arg = "signal_matrix_b")
  validate_two_signal_matrices(signal_matrix_a, signal_matrix_b)
  assert_raw_signal(signal_matrix_a, acknowledge_scaling,
                    name = "signal_matrix_a")
  assert_raw_signal(signal_matrix_b, acknowledge_scaling,
                    name = "signal_matrix_b")
  signal_matrix_a <- apply_mask_to_signal(signal_matrix_a, mask,
                                          name = "signal_matrix_a")
  signal_matrix_b <- apply_mask_to_signal(signal_matrix_b, mask,
                                          name = "signal_matrix_b")

  if (isTRUE(paired)) {
    validate_paired_matrices(signal_matrix_a, signal_matrix_b)
    diff <- signal_matrix_a - signal_matrix_b
    n    <- ncol(diff)
    mean_d <- rowMeans(diff)
    var_d  <- rowSums((diff - mean_d)^2) / (n - 1L)
    se     <- sqrt(var_d / n)
    t_vec  <- mean_d / se
    t_vec[!is.finite(t_vec)] <- 0
    return(t_vec)
  }

  n_a <- ncol(signal_matrix_a)
  n_b <- ncol(signal_matrix_b)

  mean_a <- rowMeans(signal_matrix_a)
  mean_b <- rowMeans(signal_matrix_b)

  var_a <- rowSums((signal_matrix_a - mean_a)^2) / (n_a - 1L)
  var_b <- rowSums((signal_matrix_b - mean_b)^2) / (n_b - 1L)

  se <- sqrt(var_a / n_a + var_b / n_b)
  t_vec <- (mean_a - mean_b) / se

  t_vec[!is.finite(t_vec)] <- 0
  t_vec
}

#' Validate that two signal matrices are paired (same producers)
#'
#' @keywords internal
#' @noRd
validate_paired_matrices <- function(a, b,
                                     name_a = "signal_matrix_a",
                                     name_b = "signal_matrix_b") {
  if (ncol(a) != ncol(b)) {
    cli::cli_abort(c(
      "{.arg {name_a}} and {.arg {name_b}} have different \\
       column counts.",
      "*" = "{.arg {name_a}}: {ncol(a)} producers",
      "*" = "{.arg {name_b}}: {ncol(b)} producers",
      "i" = "For a paired test, the two matrices must carry the \\
             same producers."
    ))
  }
  if (!is.null(colnames(a)) && !is.null(colnames(b))) {
    if (!identical(colnames(a), colnames(b))) {
      cli::cli_abort(c(
        "Column names of {.arg {name_a}} and {.arg {name_b}} do \\
         not match.",
        "i" = "For a paired test, column j of {.arg {name_a}} \\
               must correspond to the same producer as column j \\
               of {.arg {name_b}}."
      ))
    }
  }
  invisible(TRUE)
}
