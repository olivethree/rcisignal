#' Run every between-condition discriminability metric
#'
#' @description
#' Convenience orchestrator that runs [rel_cluster_test()] and
#' [rel_dissimilarity()] on two condition signal matrices and
#' wraps both results in an `rcisignal_rel_report` for joint
#' printing / plotting.
#'
#' Use this when you want both the spatial-pattern test and the
#' overall magnitude test in one call.
#'
#' @section Reading the result:
#' `$results$cluster_test` and `$results$dissimilarity`, one
#' result object each, fields as in the standalone functions.
#' `$method = "discriminability"`.
#'
#' @section Reliability metrics expect raw masks:
#' Both downstream metrics are scale-sensitive: the cluster test
#' uses variance-based Welch t, and Euclidean distance in
#' `rel_dissimilarity()` is sensitive to any scaling. Inputs with
#' `attr(., "source") == "rendered"` error unless
#' `acknowledge_scaling = TRUE`.
#'
#' @param signal_matrix_a,signal_matrix_b Pixels x participants,
#'   base-subtracted. Row counts must match.
#' @param img_dims Integer `c(nrow, ncol)`. If `NULL`, inferred
#'   from `attr(signal_matrix_a, "img_dims")`.
#' @param n_permutations Passed to [rel_cluster_test()]. Default
#'   2000.
#' @param n_boot Passed to [rel_dissimilarity()]. Default 2000.
#' @param cluster_threshold Passed to [rel_cluster_test()].
#'   Default 2.0.
#' @param alpha Passed to [rel_cluster_test()]. Default 0.05.
#' @param ci_level Passed to [rel_dissimilarity()]. Default 0.95.
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix_a)`. Threaded through to both downstream
#'   calls.
#' @param seed Optional integer.
#' @param progress Show `cli` progress bars.
#' @param acknowledge_scaling Logical. Forwarded.
#' @return Object of class `rcisignal_rel_report` with `$results`
#'   = named list of two result objects (`cluster_test`,
#'   `dissimilarity`) and `$method = "discriminability"`.
#' @seealso [rel_cluster_test()], [rel_dissimilarity()],
#'   [run_reliability()], [run_discriminability_pairwise()].
#' @examples
#' \dontrun{
#' # Two-condition pipeline: simulate -> per-condition CIs -> contrast.
#' sim <- simulate_briefrc_data(n_per_condition = 10, n_trials = 60,
#'                              conditions = c("target", "control"),
#'                              seed = 1)
#' a <- subset(sim$data, condition == "target")
#' b <- subset(sim$data, condition == "control")
#' sig_a <- ci_from_responses_briefrc(a, noise_matrix = sim$noise_matrix)$signal_matrix
#' sig_b <- ci_from_responses_briefrc(b, noise_matrix = sim$noise_matrix)$signal_matrix
#' run_discriminability(sig_a, sig_b,
#'                      n_permutations = 200L, n_boot = 200L, seed = 1)
#' }
#' @export
run_discriminability <- function(signal_matrix_a,
                                 signal_matrix_b,
                                 img_dims            = NULL,
                                 n_permutations      = 2000L,
                                 n_boot              = 2000L,
                                 cluster_threshold   = 2.0,
                                 alpha               = 0.05,
                                 ci_level            = 0.95,
                                 mask                = NULL,
                                 seed                = NULL,
                                 progress            = TRUE,
                                 acknowledge_scaling = FALSE) {
  validate_two_signal_matrices(signal_matrix_a, signal_matrix_b)
  if (is.null(img_dims)) {
    img_dims <- attr(signal_matrix_a, "img_dims")
  }

  results <- list(
    cluster_test = rel_cluster_test(
      signal_matrix_a, signal_matrix_b,
      img_dims            = img_dims,
      n_permutations      = n_permutations,
      cluster_threshold   = cluster_threshold,
      alpha               = alpha,
      mask                = mask,
      seed                = seed,
      progress            = progress,
      acknowledge_scaling = acknowledge_scaling
    ),
    dissimilarity = rel_dissimilarity(
      signal_matrix_a, signal_matrix_b,
      n_boot              = n_boot,
      ci_level            = ci_level,
      mask                = mask,
      seed                = seed,
      progress            = progress,
      acknowledge_scaling = acknowledge_scaling
    )
  )
  new_rcisignal_rel_report(
    results, method = "discriminability",
    img_dims = results$cluster_test$img_dims
  )
}
