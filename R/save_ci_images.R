#' Write rendered CIs to PNG or JPEG files
#'
#' @description
#' Renders each column of a `signal_matrix` (one CI) over a base face
#' image and writes the result to disk. Works for both per-producer
#' matrices (one column per producer, the `$signal_matrix` returned by
#' `ci_from_responses_*()`) and group-level matrices (one column per
#' group, the output of [group_ci()] or `$group_ci` when the generator
#' is called with `group_by =`). Filenames are derived from the column
#' names of `signal_matrix`.
#'
#' Each rendered image uses the same alpha-over compositing path as
#' [plot_ci_overlay()] (default, signed signal with blue=positive,
#' red=negative) or [plot_agreement_map()] (`palette = "fire"`,
#' unipolar `|t|`-style yellow-to-red on a per-pixel CI magnitude).
#'
#' @details
#' Filenames default to `<prefix><colname>.<ext>`, where `<prefix>` is
#' chosen automatically from `attr(signal_matrix, "ci_stage")`:
#'
#' * `"individual"` (set by `ci_from_responses_*()` on the per-producer
#'   `$signal_matrix`) -> `prefix = "ind_ci_"`.
#' * `"group"` (set by [group_ci()] on its return matrix) ->
#'   `prefix = "group_ci_"`.
#' * No `ci_stage` attribute -> defaults to `prefix = "ind_ci_"`.
#'
#' Override the auto-prefix by passing `prefix =` explicitly.
#'
#' @param signal_matrix Numeric matrix with non-empty, unique column
#'   names. Per-producer or group-level; both are accepted.
#' @param base_image Base face image. Either a numeric matrix in
#'   `[0, 1]` or a single string path to a PNG / JPEG. Used as the
#'   underlay for every rendered CI.
#' @param dir Output directory. Created (recursively) if missing.
#' @param format Output format. `"png"` (default) or `"jpeg"`.
#' @param palette Color palette. `"diverging"` (default; signed CI,
#'   blue=positive, red=negative, matching `plot_ci_overlay()`) or
#'   `"fire"` (unipolar `|t|`-style, yellow-to-red on `abs(CI)`).
#' @param prefix Optional character scalar overriding the
#'   auto-derived filename prefix (`"ind_ci_"` for per-producer
#'   matrices, `"group_ci_"` for group-level). Pass any string to
#'   force a custom convention (e.g. `prefix = "trust_"`).
#' @param threshold Optional numeric. Pixels with absolute CI value
#'   below `threshold` are forced to neutral (transparent overlay).
#' @param mask Optional logical vector of length `nrow(signal_matrix)`.
#'   Pixels with `mask = FALSE` render as base only.
#' @param zlim Optional `c(lo, hi)` color-scale endpoints. When
#'   `NULL` (default), each image is auto-scaled to its own
#'   maximum-absolute value (so colors are comparable within a
#'   producer / group, not across).
#' @param alpha_max Numeric in `[0, 1]`. Maximum opacity of the
#'   heatmap at the color-scale top. Default `0.7`.
#' @param img_dims Optional integer `c(nrow, ncol)`. Inferred from
#'   `attr(signal_matrix, "img_dims")` or from a square root of
#'   `nrow(signal_matrix)`.
#' @param quality JPEG quality in `[0, 100]`. Default `90`.
#'   Ignored for PNG.
#' @param overwrite Logical. When `FALSE` (default), the function
#'   aborts if any target file already exists. When `TRUE`, existing
#'   files are silently replaced.
#' @param quiet Logical. When `FALSE` (default), emit a one-line
#'   `cli` summary at the end.
#'
#' @return Invisibly, a character vector of the file paths written.
#'
#' @seealso [ci_from_responses_briefrc()], [ci_from_responses_2ifc()],
#'   [group_ci()], [plot_ci_overlay()], [plot_agreement_map()].
#'
#' @examples
#' \dontrun{
#' sim <- simulate_briefrc_data(
#'   n_per_condition = 10, n_trials = 60,
#'   conditions = c("A", "B"), seed = 1
#' )
#' res <- ci_from_responses_briefrc(
#'   sim$data, noise_matrix = sim$noise_matrix,
#'   base_image = sim$base_face, group_by = "condition"
#' )
#'
#' out <- tempfile("ci_export_"); dir.create(out)
#'
#' # Per-producer CIs (files: ind_ci_P001.png, ind_ci_P002.png, ...)
#' save_ci_images(res$signal_matrix, base_image = sim$base_face,
#'                dir = out)
#'
#' # Group-level CIs (files: group_ci_A.png, group_ci_B.png)
#' save_ci_images(res$group_ci, base_image = sim$base_face,
#'                dir = out, palette = "fire")
#' }
#' @export
save_ci_images <- function(signal_matrix,
                           base_image,
                           dir,
                           format    = c("png", "jpeg"),
                           palette   = c("diverging", "fire"),
                           prefix    = NULL,
                           threshold = NULL,
                           mask      = NULL,
                           zlim      = NULL,
                           alpha_max = 0.7,
                           img_dims  = NULL,
                           quality   = 90,
                           overwrite = FALSE,
                           quiet     = FALSE) {
  if (!is.matrix(signal_matrix) || !is.numeric(signal_matrix)) {
    cli::cli_abort("{.arg signal_matrix} must be a numeric matrix.")
  }
  cn <- colnames(signal_matrix)
  if (is.null(cn) || any(!nzchar(cn))) {
    cli::cli_abort(c(
      "{.arg signal_matrix} must have non-empty column names.",
      "i" = "Names become the filenames; set \\
             {.code colnames(signal_matrix)} to your producer / \\
             group ids before calling."
    ))
  }
  if (anyDuplicated(cn) > 0L) {
    cli::cli_abort(
      "{.arg signal_matrix} has duplicated column names; this \\
       would clobber output files."
    )
  }
  format  <- match.arg(format)
  palette <- match.arg(palette)
  ext <- switch(format, png = "png", jpeg = "jpeg")

  pkg <- switch(format, png = "png", jpeg = "jpeg")
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg {pkg}} is required to write {format} files.",
      "i" = "Install with {.code install.packages(\"{pkg}\")}."
    ))
  }

  if (is.null(prefix)) {
    stage <- attr(signal_matrix, "ci_stage")
    prefix <- if (identical(stage, "group")) "group_ci_" else "ind_ci_"
  }
  if (!is.character(prefix) || length(prefix) != 1L) {
    cli::cli_abort("{.arg prefix} must be a single character string.")
  }

  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (is.null(img_dims)) {
    a <- attr(signal_matrix, "img_dims")
    if (!is.null(a)) {
      img_dims <- as.integer(a)
    } else {
      side <- sqrt(nrow(signal_matrix))
      if (side != as.integer(side)) {
        cli::cli_abort(
          "Cannot infer {.arg img_dims}; pass it explicitly."
        )
      }
      img_dims <- c(as.integer(side), as.integer(side))
    }
  }
  img_dims <- as.integer(img_dims)
  if (length(img_dims) != 2L || prod(img_dims) != nrow(signal_matrix)) {
    cli::cli_abort(c(
      "{.arg img_dims} does not match {.arg signal_matrix}.",
      "*" = "img_dims: {img_dims[1L]} x {img_dims[2L]}",
      "*" = "pixels:   {nrow(signal_matrix)}"
    ))
  }

  base_mat <- resolve_base_for_overlay(base_image)
  if (!identical(as.integer(dim(base_mat)), img_dims)) {
    cli::cli_abort(c(
      "{.arg base_image} dimensions do not match {.arg signal_matrix}.",
      "*" = "base: {dim(base_mat)[1L]} x {dim(base_mat)[2L]}",
      "*" = "img_dims: {img_dims[1L]} x {img_dims[2L]}"
    ))
  }

  if (!is.null(mask)) {
    if (!is.logical(mask) || length(mask) != nrow(signal_matrix)) {
      cli::cli_abort(
        "{.arg mask} must be a logical vector of length \\
         {nrow(signal_matrix)}."
      )
    }
  }

  paths <- file.path(dir, paste0(prefix, cn, ".", ext))
  if (!isTRUE(overwrite)) {
    exists_already <- file.exists(paths)
    if (any(exists_already)) {
      n_clash <- sum(exists_already)
      cli::cli_abort(c(
        "{n_clash} output file{?s} already exist in {.path {dir}}.",
        "*" = "First: {.path {paths[exists_already][1L]}}",
        "i" = "Pass {.code overwrite = TRUE} to replace them."
      ))
    }
  }

  for (j in seq_len(ncol(signal_matrix))) {
    sig <- signal_matrix[, j]

    display <- sig
    if (!is.null(threshold)) {
      display[abs(display) < threshold] <- 0
    }
    if (!is.null(mask)) {
      display[!mask] <- 0
    }

    if (is.null(zlim)) {
      rng <- max(abs(display), na.rm = TRUE)
      if (!is.finite(rng) || rng == 0) rng <- 1
      this_zlim <- if (palette == "fire") c(0, rng) else c(-rng, rng)
    } else {
      this_zlim <- zlim
    }
    zlim_max <- max(abs(this_zlim))
    if (!is.finite(zlim_max) || zlim_max == 0) zlim_max <- 1

    norm_vec <- pmax(-1, pmin(1, display / zlim_max))
    norm_mat <- matrix(norm_vec, img_dims[1L], img_dims[2L])
    alpha_raster <- abs(norm_mat) * alpha_max

    if (palette == "fire") {
      mag_mat <- abs(norm_mat) * zlim_max
      fg_rgb  <- fire_rgb_array(mag_mat, max_mag = zlim_max)
    } else {
      fg_rgb  <- diverging_rgb_array(norm_mat)
    }
    composed <- composite_rgb_over_gray(base_mat, fg_rgb, alpha_raster)

    if (format == "png") {
      png::writePNG(composed, target = paths[j])
    } else {
      jpeg::writeJPEG(composed, target = paths[j], quality = quality / 100)
    }
  }

  if (!isTRUE(quiet)) {
    cli::cli_inform(c(
      "v" = "Wrote {ncol(signal_matrix)} {format} file{?s} to \\
             {.path {dir}}."
    ))
  }
  invisible(paths)
}
