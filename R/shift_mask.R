#' Shift a face-region mask by a number of pixels
#'
#' @description
#' Slides a logical mask by `down` and `right` pixels and returns
#' the shifted mask in the same shape as the input. Convenience
#' wrapper around the manual recipe shown in vignette section 4.5
#' for tuning an *elliptical* `make_face_mask()` region (`"full"`,
#' `"nose"`, `"mouth"`, `"upper_face"`, `"lower_face"`) on a base
#' image whose features sit a few pixels off the default geometry.
#'
#' For the *rectangle* regions (`"eyes"`, `"left_eye"`,
#' `"right_eye"`), prefer the `region_bounds` argument of
#' [make_face_mask()]: passing
#' `region_bounds = c(row_min, row_max, col_min, col_max)` in 0-1
#' image fractions is more precise than counting pixels and is
#' the canonical tuning route for those regions.
#'
#' Pixels shifted off the image are dropped, and the newly exposed
#' edge is filled with `FALSE` (i.e., outside the mask).
#'
#' @param mask Logical vector of length `prod(img_dims)`
#'   (column-major, as returned by [make_face_mask()]) or a
#'   logical matrix.
#' @param down,right Integer pixel offsets. Positive `down` moves
#'   the mask towards the bottom of the image; positive `right`
#'   moves it towards the right edge. Negative values move up /
#'   left. Defaults are `0` (no shift).
#' @param img_dims Integer `c(nrow, ncol)` (or a single integer
#'   for a square image). Required when `mask` is a vector;
#'   ignored when `mask` is already a matrix.
#' @return Logical mask of the same shape as the input
#'   (matrix-in, matrix-out; vector-in, vector-out).
#' @seealso [make_face_mask()] (the `region_bounds` route for the
#'   rectangle regions), [plot_mask_overlay()], [plot_face_mask()].
#' @export
#' @examples
#' # Reshape the default mouth mask, slide it down 20 pixels on a
#' # 256-pixel image (about 8% of the height), and confirm it
#' # still has the right number of TRUE pixels.
#' m <- make_face_mask(c(256L, 256L), region = "mouth")
#' shifted <- shift_mask(m, down = 20L, img_dims = c(256L, 256L))
#' identical(sum(m), sum(shifted))
#'
#' # Matrix-in, matrix-out.
#' mat <- matrix(m, 256L, 256L)
#' shifted_mat <- shift_mask(mat, down = 20L, right = 8L)
#' dim(shifted_mat)
shift_mask <- function(mask, down = 0L, right = 0L,
                       img_dims = NULL) {
  if (!is.numeric(down) || length(down) != 1L ||
      !is.finite(down)) {
    cli::cli_abort("{.arg down} must be a single finite number.")
  }
  if (!is.numeric(right) || length(right) != 1L ||
      !is.finite(right)) {
    cli::cli_abort("{.arg right} must be a single finite number.")
  }
  down  <- as.integer(round(down))
  right <- as.integer(round(right))

  return_vector <- FALSE
  if (is.matrix(mask)) {
    if (!is.logical(mask)) {
      cli::cli_abort("{.arg mask} matrix must be logical.")
    }
    mat <- mask
  } else if (is.logical(mask)) {
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
        "i" = "length(mask) = {length(mask)}, \\
               prod(img_dims) = {prod(img_dims)}"
      ))
    }
    mat <- matrix(mask, nrow = img_dims[1L], ncol = img_dims[2L])
    return_vector <- TRUE
  } else {
    cli::cli_abort(
      "{.arg mask} must be a logical vector or logical matrix."
    )
  }

  nr <- nrow(mat); nc <- ncol(mat)
  out <- matrix(FALSE, nr, nc)
  src_rows <- seq_len(nr) - down
  src_cols <- seq_len(nc) - right
  keep_r   <- src_rows >= 1L & src_rows <= nr
  keep_c   <- src_cols >= 1L & src_cols <= nc
  out[keep_r, keep_c] <- mat[src_rows[keep_r], src_cols[keep_c]]

  if (return_vector) as.vector(out) else out
}

#' Convert pixel-coordinate rectangle bounds to 0-1 image fractions
#'
#' @description
#' Bridge between visual inspection of a base image (where people
#' naturally think in pixel rows and columns) and the
#' `region_bounds` argument of [make_face_mask()], which expects
#' 0-1 image fractions. Returns a length-4 numeric vector ready to
#' pass straight to `make_face_mask(region_bounds = ...)`.
#'
#' @param row_min,row_max,col_min,col_max Integer pixel
#'   coordinates of the rectangle edges (1-indexed, with row 1 at
#'   the top of the image and col 1 at the left). Must satisfy
#'   `row_min <= row_max` and `col_min <= col_max`, and all four
#'   values must lie inside the image.
#' @param img_dims Integer `c(nrow, ncol)` of the target image
#'   (or a single integer for a square image).
#' @return Numeric length-4 vector
#'   `c(row_min, row_max, col_min, col_max)` in 0-1 image
#'   fractions.
#' @seealso [make_face_mask()].
#' @export
#' @examples
#' # On a 256x256 base, "rows 110-140, cols 60-115" becomes:
#' region_bounds_from_pixels(
#'   row_min = 110, row_max = 140,
#'   col_min = 60,  col_max = 115,
#'   img_dims = c(256L, 256L)
#' )
#'
#' # Pass straight to make_face_mask():
#' bounds <- region_bounds_from_pixels(110, 140, 60, 115,
#'                                     img_dims = 256L)
#' m <- make_face_mask(c(256L, 256L), region = "left_eye",
#'                     region_bounds = bounds)
#' mean(m)
region_bounds_from_pixels <- function(row_min, row_max,
                                      col_min, col_max,
                                      img_dims) {
  pix <- c(row_min, row_max, col_min, col_max)
  if (!is.numeric(pix) || length(pix) != 4L ||
      any(!is.finite(pix))) {
    cli::cli_abort(
      "{.arg row_min}, {.arg row_max}, {.arg col_min}, \\
       {.arg col_max} must each be a single finite number."
    )
  }
  img_dims <- as.integer(img_dims)
  if (length(img_dims) == 1L) img_dims <- c(img_dims, img_dims)
  if (length(img_dims) != 2L || any(img_dims < 2L)) {
    cli::cli_abort(
      "{.arg img_dims} must be an integer of length 1 or 2, \\
       each >= 2."
    )
  }
  if (row_min < 1 || row_max > img_dims[1L] ||
      col_min < 1 || col_max > img_dims[2L]) {
    cli::cli_abort(c(
      "Pixel coordinates lie outside the image.",
      "*" = "image: {img_dims[1L]} x {img_dims[2L]} (rows x cols)",
      "*" = "rows: [{row_min}, {row_max}]; cols: [{col_min}, {col_max}]"
    ))
  }
  if (row_min > row_max || col_min > col_max) {
    cli::cli_abort(
      "Require {.code row_min <= row_max} and \\
       {.code col_min <= col_max}."
    )
  }
  c(
    (row_min - 1) / (img_dims[1L] - 1),
    (row_max - 1) / (img_dims[1L] - 1),
    (col_min - 1) / (img_dims[2L] - 1),
    (col_max - 1) / (img_dims[2L] - 1)
  )
}
