#' Plot a face mask for visual verification
#'
#' Renders a face mask so you can confirm it covers the region you
#' intended before passing it to [diagnose_infoval()] or [infoval()].
#' Accepts the same input forms the diagnostics accept: a logical or
#' numeric vector (column-major, with `img_dims` supplied), a
#' logical/numeric matrix, or a path to a PNG/JPEG mask file.
#'
#' If `base_image` is supplied, the mask is drawn as a translucent
#' overlay on top of the base face — the most useful view for
#' verifying that the mask aligns with eyes, nose, mouth, etc. The
#' base image must have the same dimensions as the mask; otherwise
#' the base is dropped with a warning.
#'
#' @param mask One of: a logical or numeric vector of length
#'   `prod(img_dims)` (column-major, as returned by [make_face_mask()] or
#'   [read_face_mask()]); a logical or numeric matrix; or a path to
#'   a PNG/JPEG mask file. Pass `NULL` and a `region` to build the
#'   mask internally via [make_face_mask()].
#' @param img_dims Integer `c(nrow, ncol)`, or a single integer for a
#'   square image. Required when `mask` is a vector or when
#'   `region` is supplied; ignored otherwise.
#' @param region Optional character region name (one of the
#'   [make_face_mask()] choices) used as a convenience shortcut:
#'   builds the mask internally via
#'   `make_face_mask(img_dims, region, region_bounds)`. Mutually
#'   exclusive with `mask`.
#' @param region_bounds Optional length-4 numeric vector forwarded
#'   to [make_face_mask()] when `region` is one of the rectangle
#'   regions. Ignored otherwise.
#' @param base_image Optional path to a PNG or JPEG base face. When
#'   supplied, the mask is rendered as a translucent overlay on top.
#' @param alpha Numeric in `[0, 1]`. Overlay opacity. Default `0.5`.
#' @param col Highlight colour for the masked region. Default
#'   `"red"`.
#' @param threshold When `mask` is a numeric matrix or image path,
#'   pixels strictly above this value are treated as inside the
#'   mask. Default `0.5`. Ignored for logical input.
#' @param main Optional plot title.
#' @param ... Reserved for future use.
#'
#' @return Invisibly returns the resolved logical matrix
#'   (`nrow` x `ncol`, top-left origin).
#'
#' @seealso [plot_mask_overlay()] (overlay on a specific base image),
#'   [make_face_mask()], [read_face_mask()], [diagnose_infoval()].
#'
#' @examples
#' m <- make_face_mask(c(128L, 128L), region = "eyes")
#' plot_face_mask(m, img_dims = c(128L, 128L), main = "eyes region")
#'
#' @export
plot_face_mask <- function(mask          = NULL,
                           img_dims      = NULL,
                           base_image    = NULL,
                           alpha         = 0.5,
                           col           = "red",
                           threshold     = 0.5,
                           main          = NULL,
                           region        = NULL,
                           region_bounds = NULL,
                           ...) {
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      alpha < 0 || alpha > 1) {
    cli::cli_abort("{.arg alpha} must be a single number in [0, 1].")
  }
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
    if (is.null(img_dims)) {
      cli::cli_abort(
        "{.arg img_dims} is required when {.arg region} is supplied."
      )
    }
    mask <- make_face_mask(img_dims, region = region,
                           region_bounds = region_bounds)
  } else if (!is.null(region_bounds)) {
    cli::cli_abort(
      "{.arg region_bounds} is only valid with {.arg region}."
    )
  }

  mat <- resolve_mask_to_matrix(mask, img_dims, threshold = threshold)
  nr <- nrow(mat); nc <- ncol(mat)

  base <- NULL
  if (!is.null(base_image)) {
    if (is.character(base_image) && length(base_image) == 1L) {
      base <- read_image_as_gray(base_image)
    } else if (is.matrix(base_image) && is.numeric(base_image)) {
      base <- base_image
    } else {
      cli::cli_abort(
        "{.arg base_image} must be a file path or a numeric matrix."
      )
    }
    if (!identical(dim(base), c(nr, nc))) {
      cli::cli_warn(c(
        "Base image dims do not match mask dims; dropping base.",
        "i" = "mask: {nr}x{nc}; base: {dim(base)[1L]}x{dim(base)[2L]}"
      ))
      base <- NULL
    } else {
      base[base < 0] <- 0
      base[base > 1] <- 1
    }
  }

  op <- graphics::par(mar = c(0.5, 0.5,
                              if (is.null(main)) 0.5 else 2, 0.5))
  on.exit(graphics::par(op), add = TRUE)

  graphics::plot.new()
  graphics::plot.window(xlim = c(0, nc), ylim = c(0, nr),
                        asp = 1, xaxs = "i", yaxs = "i")

  if (!is.null(base)) {
    graphics::rasterImage(base, 0, 0, nc, nr, interpolate = FALSE)
  } else {
    graphics::rect(0, 0, nc, nr, col = "grey90", border = NA)
  }

  rgba <- grDevices::col2rgb(col) / 255
  fill <- grDevices::rgb(rgba[1L], rgba[2L], rgba[3L], alpha = alpha)
  overlay <- matrix(grDevices::rgb(0, 0, 0, 0), nr, nc)
  overlay[mat] <- fill
  graphics::rasterImage(overlay, 0, 0, nc, nr, interpolate = FALSE)

  graphics::rect(0, 0, nc, nr, border = "grey40", lwd = 1)

  if (!is.null(main)) graphics::title(main = main)

  invisible(mat)
}

# Resolve any of the accepted mask input forms (vector + img_dims,
# matrix, or path) into a logical (nr x nc) matrix with top-left
# origin — matching the orientation read_image_as_gray() and
# make_face_mask() use.
#' @keywords internal
#' @noRd
resolve_mask_to_matrix <- function(mask, img_dims, threshold = 0.5) {
  if (is.character(mask) && length(mask) == 1L) {
    img <- read_image_as_gray(mask)
    return(img > threshold)
  }
  if (is.matrix(mask)) {
    if (is.logical(mask)) return(mask)
    if (is.numeric(mask)) return(mask > threshold)
    cli::cli_abort("{.arg mask} matrix must be logical or numeric.")
  }
  if (is.logical(mask) || is.numeric(mask)) {
    if (is.null(img_dims)) {
      cli::cli_abort(
        "{.arg img_dims} is required when {.arg mask} is a vector."
      )
    }
    img_dims <- as.integer(img_dims)
    if (length(img_dims) == 1L) img_dims <- c(img_dims, img_dims)
    if (length(img_dims) != 2L || any(img_dims < 1L)) {
      cli::cli_abort(
        "{.arg img_dims} must be a positive length-1 or 2 integer."
      )
    }
    if (length(mask) != prod(img_dims)) {
      cli::cli_abort(c(
        "{.arg mask} length does not match {.arg img_dims}.",
        "i" = "length(mask) = {length(mask)}, prod(img_dims) = {prod(img_dims)}"
      ))
    }
    mat <- matrix(mask, nrow = img_dims[1L], ncol = img_dims[2L])
    if (is.numeric(mat)) mat <- mat > threshold
    return(mat)
  }
  cli::cli_abort(
    "{.arg mask} must be a logical/numeric vector or matrix, or a file path."
  )
}
