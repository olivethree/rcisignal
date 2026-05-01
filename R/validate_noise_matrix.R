#' Validate a noise matrix
#'
#' @description
#' Runs basic sanity checks on a noise matrix: numeric, finite,
#' expected dimensions. Returns an [rcisignal_diag_result()] with
#' status `"pass"`, `"warn"`, or `"fail"` rather than aborting,
#' because users typically want to see all problems at once, not
#' just the first.
#'
#' @param mat A numeric matrix, typically from [read_noise_matrix()].
#' @param expected_pixels Optional integer. If supplied, checks that
#'   `nrow(mat) == expected_pixels`. For example, 128*128 = 16384.
#' @param expected_stimuli Optional integer. If supplied, checks
#'   that `ncol(mat) == expected_stimuli`.
#' @return An object of class `"rcisignal_diag_result"`.
#' @seealso [read_noise_matrix()], [rcisignal_diag_result()].
#' @export
#' @examples
#' mat <- matrix(rnorm(16384 * 10, sd = 0.05),
#'               nrow = 16384, ncol = 10)
#' validate_noise_matrix(mat,
#'                       expected_pixels  = 16384,
#'                       expected_stimuli = 10)
validate_noise_matrix <- function(mat,
                                  expected_pixels  = NULL,
                                  expected_stimuli = NULL) {
  label  <- "Noise matrix"
  detail <- character()
  status <- "pass"

  if (!is.matrix(mat) || !is.numeric(mat)) {
    return(rcisignal_diag_result(
      "fail", label,
      "Argument is not a numeric matrix."
    ))
  }

  n_finite <- sum(is.finite(mat))
  n_total  <- length(mat)
  if (n_finite < n_total) {
    status <- "fail"
    detail <- c(
      detail,
      sprintf(
        "%d of %d entries are NA/NaN/Inf.",
        n_total - n_finite, n_total
      )
    )
  }

  if (!is.null(expected_pixels) && nrow(mat) != expected_pixels) {
    status <- "fail"
    detail <- c(
      detail,
      sprintf(
        "Expected %d rows (pixels), found %d.",
        expected_pixels, nrow(mat)
      )
    )
  }

  if (!is.null(expected_stimuli) && ncol(mat) != expected_stimuli) {
    status <- "fail"
    detail <- c(
      detail,
      sprintf(
        "Expected %d columns (stimuli), found %d.",
        expected_stimuli, ncol(mat)
      )
    )
  }

  if (status == "pass") {
    detail <- sprintf(
      "%d pixels x %d stimuli, all finite. Range: [%.3g, %.3g].",
      nrow(mat), ncol(mat), min(mat), max(mat)
    )
  }

  rcisignal_diag_result(
    status, label, detail,
    data = list(
      dim       = dim(mat),
      range     = range(mat),
      n_missing = n_total - n_finite
    )
  )
}
