#' Detect response-code inversion
#'
#' Some participants (or whole batches of exported data) may have their
#' response codes flipped relative to the convention `rcicr` expects:
#' `+1` means "inverted version chosen" when the analyst intended
#' "original chosen", or the CSV export inadvertently negated the
#' response column. A straightforward way to detect this is to compute
#' the infoVal twice per participant -- once with the original codes
#' and once with every response negated -- and compare. For correctly
#' coded data, the original should score higher; if the flipped CI
#' scores meaningfully higher, the response column is probably inverted.
#'
#' This check computes per-producer infoVal z-scores via the
#' package-native [infoval()] pipeline twice (original and
#' sign-flipped) and reports participants whose flipped infoVal
#' exceeds the original by `margin` or more. Supports both 2IFC and
#' Brief-RC. A non-zero count is a strong signal that the response
#' column is miscoded for those participants. Runs two infoVal sweeps,
#' so it takes roughly twice as long.
#'
#' @param responses Data frame with one row per trial. Required
#'   columns: `participant_id`, `stimulus`, `response` (values in
#'   `{-1, +1}`). Load yours from CSV via [read_responses()] or
#'   [utils::read.csv()]; column names are configurable via the
#'   `col_*` arguments.
#' @param method `"2ifc"` or `"briefrc"`.
#' @param rdata Path to the rcicr `.RData` file (2IFC). Either
#'   `rdata` or `stimuli` must be supplied for the 2IFC path.
#' @param stimuli In-memory stimuli list (the `$stimuli` element
#'   of an `rcisignal_sim` object). Use in place of `rdata` when
#'   the file path no longer resolves (e.g. after [saveRDS()]/
#'   [readRDS()] across R sessions).
#' @param noise_matrix Noise matrix for the Brief-RC path. Either a
#'   numeric matrix of `n_pixels x pool_size`, or a path to a text /
#'   `.rds` file (see [read_noise_matrix()]).
#' @param base_image Name of the base image in `rdata$base_face_files`
#'   (2IFC only). Default `"base"`.
#' @param col_participant,col_stimulus,col_response Column names.
#' @param margin Numeric. Flagging threshold: flipped infoVal must
#'   exceed original infoVal by at least this amount. Default `1.96`.
#' @param iter Reference-distribution iterations. Default `1000L`.
#' @param seed Optional integer; RNG state restored on exit.
#' @param ... Unused.
#'
#' @return An [rcisignal_diag_result()] object. `data$per_participant` has
#'   `participant_id`, `infoval_original`, `infoval_flipped`, and
#'   `likely_inverted` (logical).
#'
#' @examples
#' \dontrun{
#' sim <- simulate_2ifc_data(n_per_condition = 10, n_trials = 60, seed = 1)
#' check_response_inversion(sim$data, method = "2ifc",
#'                          rdata = sim$rdata_path)
#' }
#'
#' @export
check_response_inversion <- function(responses,
                                     method = c("2ifc", "briefrc"),
                                     rdata = NULL,
                                     stimuli = NULL,
                                     noise_matrix = NULL,
                                     base_image = "base",
                                     col_participant = "participant_id",
                                     col_stimulus = "stimulus",
                                     col_response = "response",
                                     margin = 1.96,
                                     iter = 1000L,
                                     seed = NULL,
                                     ...) {
  method <- match.arg(method)
  label  <- "Response inversion"

  iv_orig <- per_producer_infoval(
    responses,
    method          = method,
    rdata           = rdata,
    stimuli         = stimuli,
    noise_matrix    = noise_matrix,
    base_image      = base_image,
    col_participant = col_participant,
    col_stimulus    = col_stimulus,
    col_response    = col_response,
    iter            = iter,
    seed            = seed
  )

  flipped_responses <- responses
  flipped_responses[[col_response]] <- -flipped_responses[[col_response]]
  iv_flip <- per_producer_infoval(
    flipped_responses,
    method          = method,
    rdata           = rdata,
    stimuli         = stimuli,
    noise_matrix    = noise_matrix,
    base_image      = base_image,
    col_participant = col_participant,
    col_stimulus    = col_stimulus,
    col_response    = col_response,
    iter            = iter,
    seed            = seed
  )

  per_p <- data.frame(
    participant_id    = names(iv_orig),
    infoval_original  = unname(iv_orig),
    infoval_flipped   = unname(iv_flip[names(iv_orig)]),
    stringsAsFactors  = FALSE
  )
  per_p$delta <- per_p$infoval_flipped - per_p$infoval_original
  per_p$likely_inverted <- per_p$delta >= margin

  n_inv <- sum(per_p$likely_inverted, na.rm = TRUE)
  n_p   <- nrow(per_p)

  status <- if (n_inv == 0L) {
    "pass"
  } else if (n_inv / n_p >= 0.50) {
    "fail"
  } else {
    "warn"
  }

  detail <- c(
    sprintf(
      paste(
        "%d of %d participants have flipped infoVal exceeding",
        "original by >= %.2f."
      ),
      n_inv, n_p, margin
    ),
    if (n_inv > 0L) {
      "The flagged participants' response codes are likely inverted."
    }
  )

  rcisignal_diag_result(
    status, label, detail,
    data = list(
      per_participant = per_p,
      margin          = margin
    )
  )
}
