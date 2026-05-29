#' Run every within-condition reliability metric
#'
#' @description
#' Convenience orchestrator that runs the reliability metrics
#' proper (split-half with Spearman-Brown correction and ICC) on a
#' single condition's signal matrix and wraps the two results in
#' an `rcisignal_rel_report` for joint printing and plotting.
#'
#' Use this when you want the full within-condition reliability
#' report in one call. Pass each metric individually if you need
#' to tune arguments per metric.
#'
#' @section What is included (and what is not):
#' Returns the two metrics that quantify the reliability of the
#' group-level classification image proper: split-half (a
#' permutation-based estimate of group-CI stability with
#' Spearman-Brown projection to the full sample) and ICC(3,*) (the
#' psychometric variance decomposition). These are non-redundant.
#'
#' Leave-one-out influence screening lives in [rel_loo()] and is
#' **not** bundled here. Its output is an influence diagnostic,
#' not a reliability statistic, and mixing it into the reliability
#' report invites mis-reading `r_loo` values (which are near 1 by
#' construction) as reliability.
#'
#' @section Reading the result:
#' `$results$split_half`, `$results$icc`, one result object each,
#' with the same fields as the standalone functions.
#' `$method = "reliability"`.
#'
#' @section Reliability metrics expect raw masks:
#' Both downstream metrics expect the raw mask. Inputs with
#' `attr(., "source") == "rendered"` (set automatically by Mode 1
#' readers like [extract_signal()]) error in [rel_icc()] unless
#' `acknowledge_scaling = TRUE`.
#'
#' @param signal_matrix Pixels x participants, base-subtracted.
#' @param n_permutations Passed to [rel_split_half()]. Default 2000.
#' @param null Forwarded to [rel_split_half()]. Default `"none"`.
#' @param noise_matrix Required when
#'   `null = "random_responders"`; forwarded to [rel_split_half()].
#' @param icc_variants Passed to [rel_icc()].
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix)` (column-major) threaded through to both
#'   metrics. Build with [make_face_mask()] (parametric oval and
#'   sub-regions) or [read_face_mask()] (PNG/JPEG mask).
#' @param seed Optional integer; used for the split-half
#'   permutations.
#' @param progress Show `cli` progress bars.
#' @param acknowledge_scaling Logical. Forwarded to [rel_icc()].
#' @return Object of class `rcisignal_rel_report` with `$results`
#'   = named list of two result objects (`split_half`, `icc`)
#'   and `$method = "reliability"`.
#' @seealso [rel_split_half()], [rel_icc()], [rel_loo()] for
#'   the influence diagnostic; [run_discriminability()].
#' @examples
#' \dontrun{
#' # End-to-end demo: simulate -> compute CIs -> assess reliability.
#' sim <- simulate_briefrc_data(n_per_condition = 10, n_trials = 60, seed = 1)
#' cis <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
#' rel <- run_reliability(cis$signal_matrix,
#'                        n_permutations = 200L, seed = 1)
#' print(rel)
#' }
#' @export
run_reliability <- function(signal_matrix,
                            n_permutations      = 2000L,
                            null                = c("none",
                                                     "permutation",
                                                     "random_responders"),
                            noise_matrix        = NULL,
                            icc_variants        = c("3_1", "3_k"),
                            mask                = NULL,
                            seed                = NULL,
                            progress            = TRUE,
                            acknowledge_scaling = FALSE) {
  null <- match.arg(null)
  if (identical(null, "random_responders") && is.null(noise_matrix)) {
    cli::cli_abort(c(
      "{.arg noise_matrix} is required when \\
       {.code null = \"random_responders\"}.",
      "i" = "The null distribution simulates random-responder CIs \\
             by sign-flipping the noise pool.",
      "i" = "Pass {.arg noise_matrix} (the noise pool used at \\
             stimulus generation; load it via \\
             {.fn read_noise_matrix} if you have an rcicr rdata or \\
             a Brief-RC text file)."
    ))
  }
  validate_signal_matrix(signal_matrix)
  img_dims <- attr(signal_matrix, "img_dims")

  results <- list(
    split_half = rel_split_half(
      signal_matrix,
      n_permutations = n_permutations,
      null           = null,
      noise_matrix   = noise_matrix,
      mask           = mask,
      seed           = seed,
      progress       = progress
    ),
    icc = rel_icc(
      signal_matrix,
      variants            = icc_variants,
      mask                = mask,
      acknowledge_scaling = acknowledge_scaling
    )
  )
  new_rcisignal_rel_report(results, method = "reliability",
                          img_dims = img_dims)
}
