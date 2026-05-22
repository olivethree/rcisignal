#' Compute individual Brief-RC CIs from trial-level responses
#'
#' @description
#' Native implementation of Schmitz, Rougier & Yzerbyt (2024)'s
#' Brief-RC mask. Does not call any `rcicr::*_brief` function;
#' those do not exist in upstream rcicr v1.0.1. Only rcicr's
#' noise-pattern pool (the `stimuli` object inside an `.Rdata`
#' from `generateStimuli2IFC()`) is reused; the mask is computed
#' in pure R.
#'
#' Use this when you have Brief-RC 12 or Brief-RC 20 trial-level
#' responses and want the package to produce per-producer noise
#' masks ready for the reliability metrics.
#'
#' @details
#' Formula (Schmitz's `genMask()` exactly):
#' ```
#' X    <- data.table(response, stim)
#' X    <- X[, .(response = mean(response)), stim]  # collapse duplicates
#' mask <- (noiseMatrix[, X$stim] %*% X$response) / length(X$response)
#' ```
#' The `length(X$response)` denominator is the number of unique
#' pool ids chosen by that participant, not the raw trial count.
#' If a participant chooses the same stimulus on two trials with
#' opposite responses, those two cancel.
#'
#' The formula is symmetric in the per-trial split (6/6 for
#' Brief-RC 12, 10/10 for Brief-RC 20), so the same code path
#' handles both variants. The `method` argument is recorded as
#' provenance metadata and validated; it does not branch the
#' computation.
#'
#' @section Reading the result:
#' * `$signal_matrix` is the raw mask per producer; pass this and
#'   only this to reliability metrics or any external infoVal
#'   computation.
#' * `$rendered_ci`, when present, is `base + scaling(mask)` per
#'   producer. **Visualisation only.**
#' * `$participants` and `$img_dims` are convenience metadata.
#'
#' @section Common mistakes:
#' * Passing the "expanded" multi-row-per-trial response format.
#'   Brief-RC data is **one row per trial**; `stimulus` is the
#'   chosen pool id and `response` is `+1` (original) or `-1`
#'   (inverted). See Schmitz et al. (2024) sec. 3.1.2.
#' * Using `$rendered_ci` for downstream stats; it exists only
#'   because the user often wants to save a PNG.
#'
#' @param responses Data frame with one row per trial. Must contain
#'   the columns named by `participant_col`, `stimulus_col`,
#'   `response_col`. `response` values must be in `{-1, +1}`.
#' @param rdata_path,noise_matrix Exactly one must be supplied.
#'   Provide `rdata_path` to read the noise matrix from an rcicr
#'   `.Rdata`, or pass a pre-loaded `noise_matrix` directly.
#' @param base_image Base face image. Either a numeric matrix in
#'   `[0, 1]` (e.g. `sim$base_face` from
#'   [simulate_briefrc_data()], or a hand-built mask) or a single
#'   string path to a PNG / JPEG. Optional when `scaling = "none"`
#'   (default), in which case the image is not needed for the
#'   mathematical mask. Required for the visualisation-only
#'   `$rendered_ci` field when `scaling` is `"matched"` or
#'   `"constant"`.
#' @param base_image_path **Deprecated.** Use `base_image` (which
#'   accepts both a numeric matrix and a path). The old name still
#'   works for one release with a deprecation warning.
#' @param participant_col,stimulus_col,response_col Column names.
#' @param method Brief-RC variant: `"briefrc12"` (12 alternatives
#'   per trial, 6 original + 6 inverted; the default) or
#'   `"briefrc20"` (20 alternatives per trial, 10 original + 10
#'   inverted). Both variants are validated in Schmitz, Rougier,
#'   & Yzerbyt (2024). The CI computation is identical for both;
#'   the argument is kept as metadata so that downstream code and
#'   reports can record which paradigm produced the data.
#' @param scaling Visualisation-only scaling for the optional
#'   `$rendered_ci` field. One of `"none"` (default), `"matched"`
#'   (stretch mask to base range, then add to base) or
#'   `"constant"` (multiply mask by `scaling_constant`, then add
#'   to base). The mathematical `$signal_matrix` is always the raw
#'   unscaled mask.
#' @param scaling_constant Numeric multiplier used when
#'   `scaling = "constant"`. Ignored otherwise.
#' @return A list with `signal_matrix`, optionally `rendered_ci`,
#'   `participants`, `img_dims`, `scaling`, and `method` (the
#'   Brief-RC variant the call was made with).
#' @seealso [ci_from_responses_2ifc()], [run_reliability()],
#'   [run_discriminability()], [group_ci()] (stage 2, optional)
#' @references
#' Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the
#' brief reverse correlation: an improved tool to assess visual
#' representations. *European Journal of Social Psychology*.
#' \doi{10.1002/ejsp.3100}
#' @export
#' @examples
#' \dontrun{
#' # In your study, you replace the simulated inputs with your own:
#' #   - `sim$data` -> your responses data frame (CSV via read_responses())
#' #   - `sim$noise_matrix` -> a numeric matrix `n_pixels x pool_size`
#' #     (load via read_noise_matrix() from the OSF txt that came with
#' #     your Brief-RC stimuli).
#' # The `simulate_briefrc_data()` helper just produces the same shapes so
#' # the example below is runnable end-to-end.
#' sim <- simulate_briefrc_data(n_per_condition = 10, n_trials = 60, seed = 1)
#' res <- ci_from_responses_briefrc(
#'   sim$data,
#'   noise_matrix = sim$noise_matrix,
#'   base_image   = sim$base_face
#' )
#' dim(res$signal_matrix)   # n_pixels x n_producers
#' rel_split_half(res$signal_matrix, n_permutations = 200L, seed = 1)
#' }
ci_from_responses_briefrc <- function(responses,
                                      rdata_path       = NULL,
                                      noise_matrix     = NULL,
                                      base_image       = NULL,
                                      participant_col  = "participant_id",
                                      stimulus_col     = "stimulus",
                                      response_col     = "response",
                                      method           = c("briefrc12",
                                                           "briefrc20"),
                                      scaling          = c("none",
                                                           "matched",
                                                           "constant"),
                                      scaling_constant = NULL,
                                      base_image_path  = NULL) {
  method  <- match.arg(method)
  scaling <- match.arg(scaling)
  if (!is.null(base_image_path)) {
    cli::cli_warn(c(
      "{.arg base_image_path} is deprecated; use {.arg base_image}.",
      "i" = "{.arg base_image} accepts either a numeric matrix in \\
             {.code [0, 1]} or a path to a PNG / JPEG."
    ))
    if (is.null(base_image)) {
      base_image <- base_image_path
    }
  }
  if (scaling == "constant") {
    if (is.null(scaling_constant) ||
          !is.numeric(scaling_constant) ||
          length(scaling_constant) != 1L ||
          !is.finite(scaling_constant)) {
      cli::cli_abort(
        "{.arg scaling_constant} must be a finite numeric scalar \\
         when {.code scaling = \"constant\"}."
      )
    }
  }

  responses <- as.data.frame(responses)
  required <- c(participant_col, stimulus_col, response_col)
  missing_cols <- setdiff(required, colnames(responses))
  if (length(missing_cols) > 0L) {
    n_missing <- length(missing_cols)
    cli::cli_abort(c(
      "Missing {n_missing} column{?s} in {.arg responses}:",
      "*" = "{.val {missing_cols}}"
    ))
  }

  resp_values <- as.numeric(responses[[response_col]])
  if (!all(is.finite(resp_values))) {
    cli::cli_abort(
      "Column {.var {response_col}} contains non-finite values."
    )
  }
  uniq <- sort(unique(resp_values))
  if (!identical(uniq, c(-1, 1))) {
    msg <- c(
      "Column {.var {response_col}} must contain only values in \\
       {.val {c(-1, 1)}}.",
      "*" = "Got: {.val {uniq}}"
    )
    if (identical(uniq, c(0, 1)) ||
          identical(uniq, c(0, 1L))) {
      msg <- c(
        msg,
        "i" = "Did you mean {.code -1 / +1}? The {.code {{0, 1}}} \\
               coding is the most common silent failure in RC \\
               pipelines.",
        "i" = "Recode in one line: \\
               {.code responses${response_col} <- \\
                       2 * responses${response_col} - 1}"
      )
    }
    cli::cli_abort(msg)
  }

  if (is.null(noise_matrix) && is.null(rdata_path)) {
    cli::cli_abort(
      "Pass either {.arg rdata_path} or {.arg noise_matrix}."
    )
  }
  if (!is.null(rdata_path) && !is.null(noise_matrix)) {
    cli::cli_warn(
      "Both {.arg rdata_path} and {.arg noise_matrix} supplied; \\
       using {.arg noise_matrix}."
    )
  }
  if (is.null(noise_matrix)) {
    noise_matrix <- read_noise_matrix(rdata_path)
  } else {
    if (!is.matrix(noise_matrix)) {
      noise_matrix <- as.matrix(noise_matrix)
    }
    storage.mode(noise_matrix) <- "double"
  }

  base_resolved <- resolve_base_image(base_image, nrow(noise_matrix),
                                     scaling)
  img_dims <- base_resolved$img_dims
  base_vec <- base_resolved$base_vec
  n_pool   <- ncol(noise_matrix)

  participants <- unique(as.character(responses[[participant_col]]))
  signal_matrix <- matrix(
    NA_real_,
    nrow = nrow(noise_matrix),
    ncol = length(participants),
    dimnames = list(NULL, participants)
  )

  response_char <- as.character(responses[[participant_col]])
  stim_all      <- as.integer(responses[[stimulus_col]])
  resp_all      <- resp_values

  if (any(stim_all < 1L | stim_all > n_pool)) {
    bad <- range(stim_all)
    cli::cli_abort(c(
      "Column {.var {stimulus_col}} has ids outside the pool range.",
      "*" = "Range in data: [{bad[1]}, {bad[2]}]",
      "*" = "Noise matrix pool size: {n_pool}"
    ))
  }

  for (pid in participants) {
    idx <- which(response_char == pid)
    x <- data.table::data.table(
      response = resp_all[idx],
      stim     = stim_all[idx]
    )
    x <- x[, list(response = mean(response)), by = "stim"]
    mask <- (noise_matrix[, x$stim, drop = FALSE] %*% x$response) /
      nrow(x)
    signal_matrix[, pid] <- as.numeric(mask)
  }
  attr(signal_matrix, "img_dims") <- img_dims
  attr(signal_matrix, "source")   <- "raw"

  out <- list(
    signal_matrix = signal_matrix,
    participants  = participants,
    img_dims      = img_dims,
    scaling       = scaling,
    method        = method
  )

  if (scaling != "none") {
    out$rendered_ci <- render_brief_ci(
      signal_matrix    = signal_matrix,
      base_vec         = base_vec,
      scaling          = scaling,
      scaling_constant = scaling_constant
    )
    attr(out$rendered_ci, "img_dims") <- img_dims
    attr(out$rendered_ci, "source")   <- "rendered"
  }

  out
}

#' Resolve the `base_image` argument to a usable matrix and dims
#'
#' Accepts either a numeric matrix in `[0, 1]` (used directly) or a
#' single string path (read via the internal `read_image_as_gray()`).
#' Optional
#' (returns `base_vec = NULL`) when `scaling = "none"`, since the
#' base face only feeds the visualisation-only `$rendered_ci`. Image
#' dimensions are taken from the matrix / file when supplied; when
#' the base is omitted, dims are inferred as a square grid from
#' `nrow(noise_matrix)` and the function aborts if the pixel count
#' is not a perfect square.
#'
#' @keywords internal
#' @noRd
resolve_base_image <- function(base_image, n_pixels, scaling) {
  if (is.null(base_image)) {
    if (scaling != "none") {
      cli::cli_abort(c(
        "{.arg base_image} is required when {.code scaling != \"none\"}.",
        "i" = "Pass either a numeric matrix in {.code [0, 1]} (e.g. \\
               {.code sim$base_face}) or a file path (PNG / JPEG)."
      ))
    }
    side <- sqrt(n_pixels)
    if (side != as.integer(side)) {
      cli::cli_abort(c(
        "Cannot infer image dimensions from the noise matrix.",
        "*" = "Noise matrix has {n_pixels} pixel rows; not a square \\
               image side length.",
        "i" = "Pass {.arg base_image} so the image dims can be \\
               read directly."
      ))
    }
    side <- as.integer(side)
    return(list(img_dims = c(side, side), base_vec = NULL))
  }
  if (is.matrix(base_image) && is.numeric(base_image)) {
    if (any(!is.finite(base_image) & !is.na(base_image))) {
      cli::cli_abort("{.arg base_image} matrix contains non-finite values.")
    }
    rng_ok <- all(base_image >= 0 & base_image <= 1, na.rm = TRUE)
    if (!rng_ok) {
      cli::cli_abort(c(
        "{.arg base_image} matrix values must lie in {.code [0, 1]}.",
        "*" = "Range observed: \\
               [{min(base_image, na.rm = TRUE)}, \\
               {max(base_image, na.rm = TRUE)}]"
      ))
    }
    img_dims <- as.integer(dim(base_image))
    if (prod(img_dims) != n_pixels) {
      cli::cli_abort(c(
        "{.arg base_image} and {.arg noise_matrix} pixel counts disagree.",
        "*" = "Base: {img_dims[1]} x {img_dims[2]} = {prod(img_dims)}",
        "*" = "Noise matrix: {n_pixels} rows"
      ))
    }
    return(list(img_dims = img_dims,
                base_vec = as.vector(base_image)))
  }
  if (is.character(base_image) && length(base_image) == 1L) {
    base_img <- read_image_as_gray(base_image)
    img_dims <- as.integer(dim(base_img))
    if (prod(img_dims) != n_pixels) {
      cli::cli_abort(c(
        "{.arg base_image} and {.arg noise_matrix} pixel counts disagree.",
        "*" = "Base: {img_dims[1]} x {img_dims[2]} = {prod(img_dims)}",
        "*" = "Noise matrix: {n_pixels} rows"
      ))
    }
    return(list(img_dims = img_dims,
                base_vec = as.vector(base_img)))
  }
  cli::cli_abort(c(
    "{.arg base_image} must be one of:",
    "*" = "a numeric matrix in {.code [0, 1]}, or",
    "*" = "a single string path to a PNG / JPEG."
  ))
}

#' Render a raw mask matrix into base+scaling(mask) for visualisation
#'
#' @keywords internal
#' @noRd
render_brief_ci <- function(signal_matrix, base_vec,
                            scaling, scaling_constant = NULL) {
  base_rng <- diff(range(base_vec, finite = TRUE))
  out <- signal_matrix
  for (j in seq_len(ncol(signal_matrix))) {
    mask <- signal_matrix[, j]
    scaled <- switch(
      scaling,
      matched  = {
        mr <- diff(range(mask, finite = TRUE))
        if (mr == 0) rep(0, length(mask)) else
          mask * (base_rng / mr)
      },
      constant = mask * scaling_constant,
      mask
    )
    out[, j] <- base_vec + scaled
  }
  out
}
