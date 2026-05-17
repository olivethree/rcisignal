#' Sanity-check a face mask against a base image
#'
#' @description
#' Renders a base face with a translucent mask overlay so you can
#' verify the mask aligns with the anatomy you intended before
#' feeding it into [infoval()], [run_reliability()], or any other
#' downstream analysis. Avoids the usual hand-rolled
#' `image() + adjustcolor()` workflow.
#'
#' Companion to [plot_face_mask()] (mask alone, optional base) and
#' [plot_ci_overlay()] (signal + base): use `plot_mask_overlay()`
#' when the question is "does the mask cover the right region of
#' this specific base image?".
#'
#' @param base_image Either a numeric matrix in `[0, 1]` (grayscale,
#'   `nrow x ncol`) or a single character path to a PNG/JPEG file.
#'   PNG/JPEG paths are read via the same grayscale conversion the
#'   rest of the package uses (BT.709 luminance for RGB inputs).
#' @param mask One of: a logical or numeric vector of length
#'   `prod(dim(base_image))` (column-major, matching the convention
#'   used by [make_face_mask()] / [read_face_mask()]); a logical or
#'   numeric matrix of the same dimensions as `base_image`; or a
#'   single character path to a PNG/JPEG mask file (resolved via
#'   [read_face_mask()]). Numeric inputs are thresholded at `0.5`.
#'   Pass `NULL` (the default) and a `region` to build the mask
#'   internally via [make_face_mask()].
#' @param region Optional character region name for the
#'   convenience shortcut: passing e.g. `region = "left_eye"`
#'   builds the mask via
#'   `make_face_mask(dim(base_image), region, region_bounds)` and
#'   skips the separate construction step. Mutually exclusive with
#'   `mask`.
#' @param region_bounds Optional length-4 numeric vector forwarded
#'   to [make_face_mask()] when `region` is one of the rectangle
#'   regions (`"eyes"`, `"left_eye"`, `"right_eye"`). Ignored
#'   otherwise.
#' @param alpha Numeric in `[0, 1]`. Opacity of the mask overlay.
#'   Default `0.35`.
#' @param overlay_col Single colour for the mask overlay. Default
#'   `"red"`.
#' @param main Optional plot title (character of length 1). Default
#'   `NULL` (no title).
#' @return Invisibly `NULL`. Called for the side-effect of drawing
#'   on the active graphics device.
#' @seealso [plot_face_mask()], [plot_ci_overlay()],
#'   [make_face_mask()], [read_face_mask()].
#' @examples
#' \dontrun{
#' # Parametric oval mask plus a flat synthetic base.
#' n_side <- 128L
#' base   <- matrix(0.5, n_side, n_side)
#' eyes   <- make_face_mask(c(n_side, n_side), region = "eyes")
#' plot_mask_overlay(base, eyes)
#' }
#'
#' \dontrun{
#' # Use the simulator's base face and a sub-region to verify
#' # that the parametric mask lands on the intended anatomy.
#' sim   <- simulate_briefrc_data(
#'   n_per_condition = 5, n_trials = 10, conditions = "target",
#'   seed = 1
#' )
#' mouth <- make_face_mask(dim(sim$base_face), region = "mouth")
#' plot_mask_overlay(sim$base_face, mouth, main = "Mouth region")
#' }
#'
#' \dontrun{
#' # Mask loaded from a PNG file; pass both as paths.
#' plot_mask_overlay("path/to/base.png", "path/to/mask.png")
#' }
#'
#' # Region shortcut: skip the make_face_mask() call.
#' base <- matrix(0.5, 128L, 128L)
#' plot_mask_overlay(base, region = "left_eye")
#' @export
plot_mask_overlay <- function(base_image,
                              mask          = NULL,
                              region        = NULL,
                              region_bounds = NULL,
                              alpha         = 0.35,
                              overlay_col   = "red",
                              main          = NULL) {
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    cli::cli_abort("{.arg alpha} must be a single number in [0, 1].")
  }
  if (!is.character(overlay_col) || length(overlay_col) != 1L) {
    cli::cli_abort("{.arg overlay_col} must be a single colour name.")
  }
  if (!is.null(main) &&
      (!is.character(main) || length(main) != 1L)) {
    cli::cli_abort(
      "{.arg main} must be NULL or a single character string."
    )
  }

  base <- resolve_base_for_overlay(base_image)
  nr   <- nrow(base)
  nc   <- ncol(base)
  base[base < 0] <- 0
  base[base > 1] <- 1

  if (!is.null(region) && !is.null(mask)) {
    cli::cli_abort(
      "Pass either {.arg mask} or {.arg region}, not both."
    )
  }
  if (is.null(region) && is.null(mask)) {
    cli::cli_abort(
      "Supply one of {.arg mask} or {.arg region}."
    )
  }
  if (!is.null(region)) {
    mask <- make_face_mask(c(nr, nc), region = region,
                           region_bounds = region_bounds)
  } else if (!is.null(region_bounds)) {
    cli::cli_abort(
      "{.arg region_bounds} is only valid with {.arg region}."
    )
  }

  mask_mat <- resolve_mask_against_base(mask, nr, nc)

  op <- graphics::par(
    mar = c(0.5, 0.5, if (is.null(main)) 0.5 else 2, 0.5),
    pty = "s"
  )
  on.exit(graphics::par(op), add = TRUE)

  graphics::plot.new()
  graphics::plot.window(xlim = c(0, nc), ylim = c(0, nr),
                        asp = 1, xaxs = "i", yaxs = "i")
  graphics::rasterImage(base, 0, 0, nc, nr, interpolate = FALSE)

  rgba <- grDevices::col2rgb(overlay_col) / 255
  fill <- grDevices::rgb(rgba[1L], rgba[2L], rgba[3L], alpha = alpha)
  overlay <- matrix(grDevices::rgb(0, 0, 0, 0), nr, nc)
  overlay[mask_mat] <- fill
  graphics::rasterImage(overlay, 0, 0, nc, nr, interpolate = FALSE)

  graphics::rect(0, 0, nc, nr, border = "grey40", lwd = 1)
  if (!is.null(main)) {
    graphics::title(main = main, line = 0.5, cex.main = 1.0,
                    font.main = 1)
  }

  invisible(NULL)
}

#' @keywords internal
#' @noRd
resolve_mask_against_base <- function(mask, nr, nc) {
  n_pix <- nr * nc
  if (is.character(mask) && length(mask) == 1L) {
    vec <- read_face_mask(mask, expected_dims = c(nr, nc))
    return(matrix(vec, nrow = nr, ncol = nc))
  }
  if (is.matrix(mask)) {
    if (!identical(as.integer(dim(mask)), c(nr, nc))) {
      cli::cli_abort(c(
        "{.arg mask} dimensions do not match {.arg base_image}.",
        "*" = "base_image: {nr} x {nc}",
        "*" = "mask:       {dim(mask)[1L]} x {dim(mask)[2L]}"
      ))
    }
    if (is.logical(mask)) return(mask)
    if (is.numeric(mask)) return(mask > 0.5)
    cli::cli_abort("{.arg mask} matrix must be logical or numeric.")
  }
  if (is.logical(mask) || is.numeric(mask)) {
    if (length(mask) != n_pix) {
      cli::cli_abort(c(
        "{.arg mask} length does not match {.arg base_image}.",
        "*" = "length(mask):           {length(mask)}",
        "*" = "prod(dim(base_image)):  {n_pix}"
      ))
    }
    vec <- if (is.logical(mask)) mask else mask > 0.5
    return(matrix(vec, nrow = nr, ncol = nc))
  }
  cli::cli_abort(
    "{.arg mask} must be a logical/numeric vector or matrix, or a file path."
  )
}
