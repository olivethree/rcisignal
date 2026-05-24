## Internal helpers: validators, the rcicr soft-dep attacher,
## RNG plumbing, progress bars, session-warning state, image I/O.
## None of these are exported.

#' Attach a package without using `library()`
#'
#' `library()` inside package code triggers R CMD check's
#' "Dependence on R services not declared". We need `foreach`,
#' `tibble`, and `dplyr` **attached** (not just loaded) when
#' `rcicr::generateStimuli2IFC()` or `computeInfoVal2IFC()` run
#' because they use `%dopar%`, `tribble()`, `%>%` and `filter()` at
#' evaluation time without namespace prefixes.
#'
#' @param pkgs Character vector of package names to attach.
#' @return Invisibly `TRUE` if all attached; errors if any missing.
#' @keywords internal
#' @noRd
ensure_attached <- function(pkgs) {
  for (pkg in pkgs) {
    if (paste0("package:", pkg) %in% search()) next
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cli::cli_abort(c(
        "Package {.pkg {pkg}} is required for this operation.",
        "i" = "Install it with {.code install.packages(\"{pkg}\")}."
      ))
    }
    suppressPackageStartupMessages(attachNamespace(pkg))
  }
  invisible(TRUE)
}

#' Validate that a file path exists
#'
#' @param path Character scalar.
#' @param arg_name Argument name in the calling function, for the
#'   error message.
#' @return Invisibly returns `path` if it exists. Aborts otherwise.
#' @keywords internal
#' @noRd
validate_path <- function(path, arg_name = "path") {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    cli::cli_abort(
      "{.arg {arg_name}} must be a single, non-NA character string."
    )
  }
  if (!file.exists(path)) {
    cli::cli_abort(c(
      "File not found: {.path {path}}",
      "i" = "Passed as {.arg {arg_name}}."
    ))
  }
  invisible(path)
}

#' Validate that a responses data frame has the expected columns
#'
#' @param responses A `data.frame` or `data.table`.
#' @param col_participant,col_stimulus,col_response Required column
#'   names.
#' @param col_rt Optional column name. If non-NULL, must be present.
#' @return Invisibly returns `responses`. Aborts if validation fails.
#' @keywords internal
#' @noRd
validate_responses_df <- function(responses,
                                  col_participant,
                                  col_stimulus,
                                  col_response,
                                  col_rt = NULL) {
  if (!is.data.frame(responses)) {
    cli::cli_abort("{.arg responses} must be a data frame.")
  }
  required <- c(col_participant, col_stimulus, col_response)
  if (!is.null(col_rt)) required <- c(required, col_rt)
  missing_cols <- setdiff(required, names(responses))
  if (length(missing_cols) > 0L) {
    n_missing <- length(missing_cols)
    cli::cli_abort(c(
      "{.arg responses} is missing {n_missing} required column{?s}: \\
       {.val {missing_cols}}",
      "i" = "Available columns: {.val {names(responses)}}"
    ))
  }
  invisible(responses)
}

#' Assert a pixels x participants signal matrix is well-formed
#'
#' Every `rel_*` function takes a signal matrix. This centralises the
#' common validation: numeric matrix, at least 4 participants (minimum
#' meaningful for split-half), warn at < 30.
#'
#' @param x Candidate signal matrix.
#' @param name Argument name, for diagnostic messages.
#' @param min_participants Hard floor below which we abort.
#' @return Invisibly the input.
#' @keywords internal
#' @noRd
validate_signal_matrix <- function(x,
                                   name = "signal_matrix",
                                   min_participants = 4L) {
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort(c(
      "{.arg {name}} must be a numeric matrix.",
      "i" = "Got {.cls {class(x)}}."
    ))
  }
  if (ncol(x) < min_participants) {
    cli::cli_abort(c(
      "{.arg {name}} must have at least {min_participants} participants \\
       (columns).",
      "i" = "Got {ncol(x)}. Reliability metrics are undefined below this."
    ))
  }
  if (ncol(x) < 30L) {
    cli::cli_warn(c(
      "{.arg {name}} has fewer than 30 participants ({ncol(x)}).",
      "i" = "Reliability estimates will be noisy; a sample size of \\
             N >= 60 is recommended for stable assessment."
    ))
  }
  invisible(x)
}

#' Validate and apply a region mask to a signal matrix
#'
#' Centralised handling for the `mask` argument now accepted by every
#' `rel_*()` function. Validates that `mask` is a logical vector of
#' the right length (matching `nrow(signal_matrix)`), then row-
#' subsets the signal matrix to the masked pixels.
#'
#' Returns the (possibly-subsetted) matrix, with the
#' `img_dims` attribute preserved when no subsetting happens (when
#' `mask = NULL`) and dropped when the mask reshapes the rows
#' (subsetted matrices are no longer 2D-image-shaped). The
#' `source` attribute (`"raw"` / `"rendered"`) is preserved across
#' subsetting because masking pixels does not change whether the
#' values are raw or scaled.
#'
#' @param signal_matrix Pixels x participants numeric matrix.
#' @param mask Optional logical vector of length `nrow(signal_matrix)`.
#' @param name Argument name of the matrix at the call site, for
#'   error messages.
#' @return The signal matrix unchanged (`mask = NULL`) or
#'   `signal_matrix[mask, , drop = FALSE]`.
#' @keywords internal
#' @noRd
apply_mask_to_signal <- function(signal_matrix, mask = NULL,
                                 name = "signal_matrix") {
  if (is.null(mask)) return(signal_matrix)
  if (!is.logical(mask)) {
    cli::cli_abort(
      "{.arg mask} must be a logical vector (got {.cls {class(mask)}})."
    )
  }
  if (length(mask) != nrow(signal_matrix)) {
    cli::cli_abort(c(
      "{.arg mask} length must match {.code nrow({name})}.",
      "*" = "{.arg mask}: {length(mask)}",
      "*" = "{.code nrow({name})}: {nrow(signal_matrix)}"
    ))
  }
  if (sum(mask) < 4L) {
    cli::cli_abort(c(
      "{.arg mask} selects too few pixels ({sum(mask)}).",
      "i" = "Reliability statistics on a few-pixel subset are \\
             dominated by noise. Use a more permissive mask."
    ))
  }
  src <- attr(signal_matrix, "source")
  out <- signal_matrix[mask, , drop = FALSE]
  attr(out, "img_dims") <- NULL
  if (!is.null(src)) attr(out, "source") <- src
  out
}

#' Assert that a signal matrix is raw (not rendered/scaled)
#'
#' The shared enforcement point for variance-based metrics
#' (`pixel_t_test()`, `rel_icc()`, `rel_dissimilarity()`).
#' `rel_cluster_test()` inherits enforcement via its internal
#' `pixel_t_test()` call.
#'
#' Decision tree:
#' * `attr(x, "source") == "raw"`: pass.
#' * `attr(x, "source") == "rendered"`: error unless
#'   `acknowledge_scaling = TRUE`.
#' * No `source` attribute: fall back to `looks_scaled()` heuristic;
#'   if flagged, emit the once-per-session warning (do not error).
#'
#' Why not a single enforcement point in `pixel_t_test()`:
#' `rel_icc()` and `rel_dissimilarity()` do not call
#' `pixel_t_test()`, so a single-point design would silently fail
#' to enforce on those code paths. Each variance-based metric calls
#' this helper directly.
#'
#' @keywords internal
#' @noRd
assert_raw_signal <- function(x,
                              acknowledge_scaling = FALSE,
                              name = "signal_matrix") {
  src <- attr(x, "source")
  if (identical(src, "raw")) return(invisible(TRUE))
  if (identical(src, "rendered") && !isTRUE(acknowledge_scaling)) {
    cli::cli_abort(c(
      "{.arg {name}} is a {.strong rendered} CI (PNG-derived); \\
       variance-based metrics give wrong numbers on rendered data.",
      "i" = "Mode 2 ({.fn ci_from_responses_2ifc} / \\
             {.fn ci_from_responses_briefrc}) returns the raw mask.",
      "i" = "If you understand the trade-off and want to proceed \\
             anyway, pass {.code acknowledge_scaling = TRUE}.",
      "i" = "See the rcisignal output-reliability vignette for the \\
             full discussion."
    ))
  }
  if (identical(src, "rendered") && isTRUE(acknowledge_scaling)) {
    return(invisible(TRUE))
  }
  if (!isTRUE(acknowledge_scaling) && looks_scaled(x)) {
    warn_looks_scaled(name)
  }
  invisible(TRUE)
}

#' Assert two signal matrices have compatible shape for between-condition
#'
#' @keywords internal
#' @noRd
validate_two_signal_matrices <- function(a, b,
                                         name_a = "signal_matrix_a",
                                         name_b = "signal_matrix_b") {
  validate_signal_matrix(a, name_a)
  validate_signal_matrix(b, name_b)
  if (nrow(a) != nrow(b)) {
    cli::cli_abort(c(
      "Pixel counts differ between conditions.",
      "*" = "{.arg {name_a}}: {nrow(a)} pixels",
      "*" = "{.arg {name_b}}: {nrow(b)} pixels"
    ))
  }
  invisible(TRUE)
}

#' Validate and normalise an `img_dims` argument
#'
#' Accepts `c(nrow, ncol)` or a single integer (interpreted as
#' square `c(n, n)`). Returns integer length-2.
#'
#' @keywords internal
#' @noRd
validate_img_dims <- function(img_dims, n_pixels,
                              name = "img_dims") {
  if (length(img_dims) == 1L) {
    img_dims <- c(img_dims, img_dims)
  }
  if (length(img_dims) != 2L || any(!is.finite(img_dims)) ||
        any(img_dims < 1L)) {
    cli::cli_abort(
      "{.arg {name}} must be a positive integer of length 1 or 2 \\
       (`c(nrow, ncol)`)."
    )
  }
  img_dims <- as.integer(img_dims)
  if (prod(img_dims) != n_pixels) {
    cli::cli_abort(c(
      "{.arg {name}} is inconsistent with the signal matrix.",
      "i" = "{img_dims[1]} x {img_dims[2]} = {prod(img_dims)} != \\
             {n_pixels} pixels."
    ))
  }
  img_dims
}

#' Optionally call `set.seed()` inside a function without leaking state
#'
#' Saves the caller's RNG state, sets the package seed for the
#' duration of `expr`, and restores the state on exit. If `seed` is
#' `NULL`, `expr` runs with the caller's current RNG.
#'
#' @keywords internal
#' @noRd
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(
      rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
      add = TRUE
    )
  }
  set.seed(seed)
  force(expr)
}

# `cli` ties a progress bar's lifetime to the environment in which it
# was created. Pin every bar to a dedicated package-private env so
# start/tick/done agree on the bar id even when a loop runs inside a
# captured expression (e.g. with_seed({...})) or under knitr.
.rcisignal_progress_env <- new.env(parent = baseenv())

#' @keywords internal
#' @noRd
progress_start <- function(total, label, show = TRUE) {
  if (!isTRUE(show)) return(NULL)
  cli::cli_progress_bar(label, total = total, clear = TRUE,
                        .envir = .rcisignal_progress_env)
}

#' @keywords internal
#' @noRd
progress_tick <- function(id) {
  if (is.null(id)) return(invisible())
  cli::cli_progress_update(id = id, .envir = .rcisignal_progress_env)
}

#' @keywords internal
#' @noRd
progress_done <- function(id) {
  if (is.null(id)) return(invisible())
  cli::cli_progress_done(id = id, .envir = .rcisignal_progress_env)
}

## ---- session-state warning helpers ----------------------------------
##
## The "raw vs rendered" caveat is the load-bearing distinction in
## this package: PNG-derived signals contain `scaling(noise)`, not
## raw `noise`. Warnings fire once per session so the user sees them
## without being hammered every call. State lives in a package-
## private env that the test suite resets via
## `reset_session_warnings()`.

.rcisignal_session_state <- new.env(parent = emptyenv())
.rcisignal_session_state$mode1_warning_emitted          <- FALSE
.rcisignal_session_state$looks_scaled_warning_emitted   <- FALSE
.rcisignal_session_state$icc_resolution_warning_emitted <- FALSE
.rcisignal_session_state$loo_sd_deprecated_emitted      <- FALSE
.rcisignal_session_state$noise_cache_built_emitted      <- FALSE
.rcisignal_session_state$noise_cache_reused_emitted     <- FALSE

#' @keywords internal
#' @noRd
reset_session_warnings <- function() {
  .rcisignal_session_state$mode1_warning_emitted          <- FALSE
  .rcisignal_session_state$looks_scaled_warning_emitted   <- FALSE
  .rcisignal_session_state$icc_resolution_warning_emitted <- FALSE
  .rcisignal_session_state$loo_sd_deprecated_emitted      <- FALSE
  .rcisignal_session_state$noise_cache_built_emitted      <- FALSE
  .rcisignal_session_state$noise_cache_reused_emitted     <- FALSE
  invisible()
}

#' Emit the once-per-session rel_loo(flag_method = "sd") deprecation
#'
#' Default is `"mad"`; `"sd"` retained for backwards compatibility
#' and slated for removal in v0.2.0.
#'
#' @keywords internal
#' @noRd
warn_loo_sd_deprecated <- function() {
  if (isTRUE(getOption("rcisignal.silence_loo_deprecation", FALSE))) {
    return(invisible())
  }
  if (isTRUE(.rcisignal_session_state$loo_sd_deprecated_emitted)) {
    return(invisible())
  }
  cli::cli_warn(c(
    "{.code rel_loo(flag_method = \"sd\")} is deprecated and will \\
     be removed in v0.2.0.",
    "i" = "The MAD/median rule (the new default) is robust to the \\
           influential producers LOO is designed to flag.",
    "*" = "Drop the {.arg flag_method} argument or pass \\
           {.code flag_method = \"mad\"} explicitly.",
    "i" = "Silence: \\
           {.code options(rcisignal.silence_loo_deprecation = TRUE)}."
  ))
  .rcisignal_session_state$loo_sd_deprecated_emitted <- TRUE
  invisible()
}

#' Emit the once-per-session ICC(3,k) resolution-asymptote warning
#'
#' Called by `rel_icc()` when the image has more than 50,000 pixels.
#' Flags that ICC(3,k) approaches 1 at large image sizes and is not
#' resolution-comparable.
#'
#' @keywords internal
#' @noRd
warn_icc_large_image <- function(n_targets) {
  if (isTRUE(getOption("rcisignal.silence_icc_warning", FALSE))) {
    return(invisible())
  }
  if (isTRUE(.rcisignal_session_state$icc_resolution_warning_emitted)) {
    return(invisible())
  }
  cli::cli_warn(c(
    "ICC(3,k) tends toward 1 at large image sizes \\
     ({n_targets} pixels).",
    "i" = "This is a property of ICC(*,k), not a reliability \\
           statement: the error mean square scales as \\
           1/((n-1)(k-1)) while the row mean square scales as \\
           1/(n-1), so ICC(3,k) asymptotes to 1 as n grows.",
    "*" = "Report {.strong ICC(3,1)} as the primary, \\
           resolution-comparable statistic when comparing across \\
           image sizes.",
    "i" = "Silence: {.code options(rcisignal.silence_icc_warning = TRUE)}."
  ))
  .rcisignal_session_state$icc_resolution_warning_emitted <- TRUE
  invisible()
}

#' Emit the once-per-session Mode-1 scaling warning
#'
#' Called by `read_cis()`, `extract_signal()`, `read_signal_matrix()`.
#' Silent if `acknowledge_scaling = TRUE`, if the option
#' `rcisignal.silence_scaling_warning` is `TRUE`, or if the warning
#' has already fired in this session.
#'
#' @keywords internal
#' @noRd
warn_mode1_scaling <- function(acknowledge_scaling = FALSE) {
  if (isTRUE(acknowledge_scaling)) return(invisible())
  if (isTRUE(getOption("rcisignal.silence_scaling_warning", FALSE))) {
    return(invisible())
  }
  if (isTRUE(.rcisignal_session_state$mode1_warning_emitted)) {
    return(invisible())
  }
  cli::cli_warn(c(
    "PNG-derived signal matrix contains the {.strong rendered} CI, \\
     not the raw mask.",
    "i" = "PNGs encode {.code base + scaling(mask)}, so \\
           {.code cis - base} recovers {.code scaling(mask)}, not \\
           the raw mask itself.",
    "*" = "Pearson-based metrics ({.fn rel_split_half}, \\
           {.fn rel_loo}) survive a uniform linear scaling but are \\
           still distorted by per-CI {.val matched} scaling.",
    "*" = "{.fn rel_icc}, Euclidean half of {.fn rel_dissimilarity}, \\
           {.fn pixel_t_test}, {.fn rel_cluster_test} and any \\
           hand-rolled {.code infoVal} computation are sensitive to \\
           {.strong any} scaling.",
    "i" = "{.code rcicr::computeInfoVal2IFC()} extracts the raw \\
           {.field $ci} component internally, so the standard 2IFC \\
           infoVal path is unaffected.",
    "i" = "Prefer Mode 2 ({.fn ci_from_responses_2ifc} / \\
           {.fn ci_from_responses_briefrc}) when raw responses are \\
           available; those return the raw mask.",
    "i" = "Silence: pass {.code acknowledge_scaling = TRUE}, set \\
           {.code options(rcisignal.silence_scaling_warning = TRUE)}, \\
           or generate PNGs with {.code scaling = \"none\"} so the \\
           output is effectively raw."
  ))
  .rcisignal_session_state$mode1_warning_emitted <- TRUE
  invisible()
}

#' Emit the once-per-session "looks scaled" warning
#'
#' Called by every `rel_*()` and `run_*()` on entry when
#' `looks_scaled()` flags the input. Quieter than the Mode-1 warning
#' (which targets the input boundary directly).
#'
#' @keywords internal
#' @noRd
warn_looks_scaled <- function(name = "signal_matrix") {
  if (isTRUE(getOption("rcisignal.silence_scaling_warning", FALSE))) {
    return(invisible())
  }
  if (isTRUE(.rcisignal_session_state$looks_scaled_warning_emitted)) {
    return(invisible())
  }
  cli::cli_warn(c(
    "{.arg {name}} looks like it may be a {.strong rendered} CI \\
     (per-column dynamic range is highly heterogeneous).",
    "i" = "If the matrix came from {.fn read_cis} / \\
           {.fn extract_signal} on rendered PNGs, downstream metrics \\
           may be distorted. Mode 2 ({.fn ci_from_responses_2ifc} / \\
           {.fn ci_from_responses_briefrc}) returns the raw mask.",
    "i" = "Heuristic only; silence with \\
           {.code options(rcisignal.silence_scaling_warning = TRUE)} \\
           if the matrix is genuinely raw."
  ))
  .rcisignal_session_state$looks_scaled_warning_emitted <- TRUE
  invisible()
}

#' Inform the user (once per session) that an .rds cache was built or reused
#'
#' Called from auto-detecting readers such as `read_noise_matrix()`
#' when they materialise a slow-to-parse text source into a sibling
#' `.rds`. Silent caching breeds confusion when the source file
#' changes; surfacing it once per session is the right balance.
#'
#' @keywords internal
#' @noRd
inform_cache_built <- function(rds_path) {
  if (isTRUE(getOption("rcisignal.silence_cache_messages", FALSE))) {
    return(invisible())
  }
  if (isTRUE(.rcisignal_session_state$noise_cache_built_emitted)) {
    return(invisible())
  }
  cli::cli_inform(c(
    "v" = "Built fast cache: {.path {rds_path}}",
    "i" = "Subsequent calls with the same source will load from the \\
           cache automatically. Delete the {.path .rds} to force a \\
           fresh parse, or silence this message with \\
           {.code options(rcisignal.silence_cache_messages = TRUE)}."
  ))
  .rcisignal_session_state$noise_cache_built_emitted <- TRUE
  invisible()
}

#' @keywords internal
#' @noRd
inform_cache_reused <- function(rds_path) {
  if (isTRUE(getOption("rcisignal.silence_cache_messages", FALSE))) {
    return(invisible())
  }
  if (isTRUE(.rcisignal_session_state$noise_cache_reused_emitted)) {
    return(invisible())
  }
  cli::cli_inform(c(
    "i" = "Loaded cached matrix from {.path {rds_path}}.",
    "i" = "Silence: \\
           {.code options(rcisignal.silence_cache_messages = TRUE)}."
  ))
  .rcisignal_session_state$noise_cache_reused_emitted <- TRUE
  invisible()
}

#' Read and convert an image file to a grayscale numeric matrix
#'
#' Handles PNG and JPEG via the optional `png` / `jpeg` packages,
#' collapses RGB channels to luminance via ITU-R BT.709 weights if
#' needed, returns an `nrow x ncol` numeric matrix with values in
#' the 0-1 range.
#'
#' @keywords internal
#' @noRd
read_image_as_gray <- function(path) {
  if (!file.exists(path)) {
    cli::cli_abort("Image file not found: {.path {path}}")
  }
  ext <- tolower(tools::file_ext(path))
  img <- switch(
    ext,
    png = {
      if (!requireNamespace("png", quietly = TRUE)) {
        cli::cli_abort(
          "Reading PNG files requires the {.pkg png} package."
        )
      }
      png::readPNG(path)
    },
    jpg = ,
    jpeg = {
      if (!requireNamespace("jpeg", quietly = TRUE)) {
        cli::cli_abort(
          "Reading JPEG files requires the {.pkg jpeg} package."
        )
      }
      jpeg::readJPEG(path)
    },
    cli::cli_abort(
      "Unsupported image extension {.val {ext}} for {.path {path}}."
    )
  )
  if (length(dim(img)) == 2L) {
    return(img)
  }
  nch <- dim(img)[3]
  if (nch >= 3L) {
    0.2126 * img[, , 1] + 0.7152 * img[, , 2] + 0.0722 * img[, , 3]
  } else {
    img[, , 1]
  }
}

#' Resolve a base-image argument (matrix or path) to a grayscale matrix
#'
#' @keywords internal
#' @noRd
resolve_base_for_overlay <- function(base_image) {
  if (is.character(base_image) && length(base_image) == 1L) {
    return(read_image_as_gray(base_image))
  }
  if (is.matrix(base_image) && is.numeric(base_image)) {
    return(base_image)
  }
  cli::cli_abort(
    "{.arg base_image} must be a numeric matrix or a path to PNG/JPEG."
  )
}

#' Build a sign-based diverging-palette RGB array
#'
#' Positive entries map to blue, negative to red, zero to (0, 0, 0).
#' Caller multiplies the result by an alpha raster when compositing.
#' Matches the color convention of `plot_ci_overlay()` and the
#' cluster-test plots (positive = blue, negative = red).
#'
#' @keywords internal
#' @noRd
diverging_rgb_array <- function(signed_mat) {
  pos <- signed_mat > 0
  neg <- signed_mat < 0
  nr  <- nrow(signed_mat)
  nc  <- ncol(signed_mat)
  fg  <- array(0, dim = c(nr, nc, 3L))
  fg[, , 1] <- ifelse(pos, 0.10, ifelse(neg, 0.85, 0))
  fg[, , 2] <- ifelse(pos, 0.20, ifelse(neg, 0.10, 0))
  fg[, , 3] <- ifelse(pos, 0.85, ifelse(neg, 0.10, 0))
  fg
}

#' Build a YlOrRd-sampled RGB array from a non-negative magnitude matrix
#'
#' `mag_mat` must be non-negative (typically `abs(t)` after thresholding).
#' Values are scaled by `max_mag` (the color-scale top) and sampled from
#' the reversed `YlOrRd` ramp. Pixels at zero get the pale-yellow end of
#' the ramp, which is intentionally near-white so the base shows through
#' when the alpha is also low.
#'
#' @keywords internal
#' @noRd
fire_rgb_array <- function(mag_mat, max_mag) {
  if (!is.finite(max_mag) || max_mag <= 0) max_mag <- 1
  nr <- nrow(mag_mat)
  nc <- ncol(mag_mat)
  ramp <- grDevices::hcl.colors(256L, "YlOrRd", rev = TRUE)
  idx  <- floor(pmin(pmax(mag_mat / max_mag, 0), 1) * 255) + 1L
  idx[is.na(idx)] <- 1L
  rgb_mat <- grDevices::col2rgb(ramp[idx]) / 255
  fg <- array(0, dim = c(nr, nc, 3L))
  fg[, , 1] <- matrix(rgb_mat[1L, ], nrow = nr, ncol = nc)
  fg[, , 2] <- matrix(rgb_mat[2L, ], nrow = nr, ncol = nc)
  fg[, , 3] <- matrix(rgb_mat[3L, ], nrow = nr, ncol = nc)
  fg
}

#' Alpha-over composite a foreground RGB raster onto a grayscale base
#'
#' Standard over-compositing: per-pixel
#' `out = (1 - alpha) * base + alpha * fg`, replicated across the
#' three RGB channels and clamped to `[0, 1]`. Used by
#' `plot_ci_overlay()`, `plot_agreement_map()`, and the cluster-test
#' plot methods so they share a single composition path.
#'
#' @keywords internal
#' @noRd
composite_rgb_over_gray <- function(base_gray, fg_rgb_array, alpha_raster) {
  nr <- nrow(base_gray)
  nc <- ncol(base_gray)
  base_layer <- array(rep(as.vector(base_gray), 3L), dim = c(nr, nc, 3L))
  composed <- array(0, dim = c(nr, nc, 3L))
  composed[, , 1] <- (1 - alpha_raster) * base_layer[, , 1] +
    alpha_raster * fg_rgb_array[, , 1]
  composed[, , 2] <- (1 - alpha_raster) * base_layer[, , 2] +
    alpha_raster * fg_rgb_array[, , 2]
  composed[, , 3] <- (1 - alpha_raster) * base_layer[, , 3] +
    alpha_raster * fg_rgb_array[, , 3]
  composed[composed < 0] <- 0
  composed[composed > 1] <- 1
  composed
}
