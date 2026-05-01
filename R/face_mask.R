#' Build an oval face-region mask for a square image
#'
#' @description
#' Returns a logical vector of length `prod(img_dims)` marking an
#' elliptical face region (or a face sub-region) centred on the
#' image. Pass to [infoval()] to restrict both observed and reference
#' Frobenius norms to the masked region; pass to any `rel_*()`
#' function via row-subsetting (`signal_matrix[mask, ]`) to compute
#' reliability or discriminability on a single anatomical region.
#'
#' Six regions are supported:
#' * `"full"` (default): the full face oval. Defaults follow the
#'   geometry used by Schmitz, Rougier, & Yzerbyt (2024) and tunable
#'   via `half_width` / `half_height`.
#' * `"eyes"`: two small ellipses at the typical eye positions.
#' * `"nose"`: a narrow vertical ellipse along the midline.
#' * `"mouth"`: a wide-and-short ellipse below centre.
#' * `"upper_face"`, `"lower_face"`: top and bottom halves of the
#'   full face oval, useful for crude top/bottom comparisons.
#'
#' Sub-region geometries are heuristic approximations matched to a
#' typical centred face on a square base (e.g. 256x256 KDEF male).
#' For non-default base images, tune `centre`, `half_width`,
#' `half_height`; sub-region positions scale relative to the
#' full-face ellipse.
#'
#' @param img_dims Integer `c(nrow, ncol)`.
#' @param region Character. One of `"full"`, `"eyes"`, `"nose"`,
#'   `"mouth"`, `"upper_face"`, `"lower_face"`.
#' @param centre Numeric `c(row, col)` in the 0-1 coordinate range.
#'   Default `c(0.5, 0.5)`.
#' @param half_width Full-face ellipse horizontal half-axis, as a
#'   fraction of image width. Default 0.35.
#' @param half_height Full-face ellipse vertical half-axis, as a
#'   fraction of image height. Default 0.45.
#' @return Logical vector of length `prod(img_dims)`, column-major
#'   to match the package's image vectorisation convention.
#' @seealso [read_face_mask()], [infoval()], [plot_agreement_map()]
#' @references
#' Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the
#' brief reverse correlation: an improved tool to assess visual
#' representations. *European Journal of Social Psychology*.
#' \doi{10.1002/ejsp.3100}
#' @export
#' @examples
#' full  <- make_face_mask(c(128L, 128L))
#' eyes  <- make_face_mask(c(128L, 128L), region = "eyes")
#' mouth <- make_face_mask(c(128L, 128L), region = "mouth")
#' c(full = mean(full), eyes = mean(eyes), mouth = mean(mouth))
make_face_mask <- function(img_dims,
                           region      = c("full", "eyes", "nose",
                                           "mouth", "upper_face",
                                           "lower_face"),
                           centre      = c(0.5, 0.5),
                           half_width  = 0.35,
                           half_height = 0.45) {
  region   <- match.arg(region)
  img_dims <- as.integer(img_dims)
  if (length(img_dims) == 1L) img_dims <- c(img_dims, img_dims)
  if (length(img_dims) != 2L || any(img_dims < 1L)) {
    cli::cli_abort(
      "{.arg img_dims} must be a positive length-1 or 2 integer."
    )
  }

  nr <- img_dims[1L]; nc <- img_dims[2L]
  full_oval <- ellipse_mask(img_dims, centre[1L], centre[2L],
                            half_height, half_width)

  if (region == "full") {
    return(as.vector(full_oval))
  }

  result <- matrix(FALSE, nr, nc)

  if (region == "eyes") {
    eye_row    <- centre[1L] - 0.10 * (2 * half_height)
    eye_offset <- 0.40 * half_width
    eye_hh     <- 0.10 * half_height
    eye_hw     <- 0.16 * half_width
    result <- result |
      ellipse_mask(img_dims, eye_row, centre[2L] - eye_offset,
                   eye_hh, eye_hw) |
      ellipse_mask(img_dims, eye_row, centre[2L] + eye_offset,
                   eye_hh, eye_hw)
  } else if (region == "nose") {
    nose_row <- centre[1L] + 0.05 * (2 * half_height)
    nose_hh  <- 0.20 * half_height
    nose_hw  <- 0.12 * half_width
    result <- ellipse_mask(img_dims, nose_row, centre[2L],
                           nose_hh, nose_hw)
  } else if (region == "mouth") {
    mouth_row <- centre[1L] + 0.32 * (2 * half_height)
    mouth_hh  <- 0.10 * half_height
    mouth_hw  <- 0.30 * half_width
    result <- ellipse_mask(img_dims, mouth_row, centre[2L],
                           mouth_hh, mouth_hw)
  } else if (region == "upper_face") {
    rr <- (row(matrix(0, nr, nc)) - 1L) / (nr - 1L)
    result <- full_oval & (rr < centre[1L])
  } else if (region == "lower_face") {
    rr <- (row(matrix(0, nr, nc)) - 1L) / (nr - 1L)
    result <- full_oval & (rr >= centre[1L])
  }

  result <- result & full_oval
  as.vector(result)
}

#' Read an image-based face-region mask from a PNG or JPEG
#'
#' @description
#' Reads a binary mask image from disk and returns a logical vector
#' of length `prod(img_dims)` in column-major order, the format
#' the rest of the package expects. White / light pixels (above
#' `threshold`) become `TRUE`, dark pixels become `FALSE`. Use this
#' for image-based masks created with
#' [webmorphR::mask_oval()](https://debruine.github.io/webmorphR/),
#' painted in GIMP, drawn with PowerPoint shapes, or any other tool
#' that produces a binary PNG/JPEG.
#'
#' Companion to [make_face_mask()] (parametric oval). Either
#' function returns the same logical-vector shape.
#'
#' @param path Path to a PNG or JPEG mask image.
#' @param threshold Numeric in `[0, 1]`. Pixels with luminance above
#'   this become `TRUE`. Default `0.5` (mid-grey).
#' @param invert If `TRUE`, the mask is inverted (useful for
#'   black-on-white masks). Default `FALSE`.
#' @param expected_dims Optional integer `c(nrow, ncol)`. When set,
#'   aborts if the loaded image's dimensions do not match — useful
#'   for catching a wrong-resolution mask before it silently corrupts
#'   a downstream computation.
#' @return Logical vector of length `prod(img_dims)`, column-major.
#' @seealso [make_face_mask()], [infoval()].
#' @export
#' @examples
#' \dontrun{
#' fm <- read_face_mask("masks/oval_256.png",
#'                      expected_dims = c(256L, 256L))
#' iv <- infoval(signal_matrix, noise_matrix, trial_counts,
#'               mask = fm, iter = 1000L, seed = 1L)
#' }
read_face_mask <- function(path, threshold = 0.5,
                           invert = FALSE,
                           expected_dims = NULL) {
  if (!file.exists(path)) {
    cli::cli_abort("Mask file not found: {.path {path}}")
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
        !is.finite(threshold) || threshold < 0 || threshold > 1) {
    cli::cli_abort(
      "{.arg threshold} must be a finite numeric in {.code [0, 1]}."
    )
  }
  ext <- tolower(tools::file_ext(path))
  img <- switch(
    ext,
    png = {
      if (!requireNamespace("png", quietly = TRUE)) {
        cli::cli_abort(
          "Reading PNG masks requires the {.pkg png} package."
        )
      }
      png::readPNG(path)
    },
    jpg = ,
    jpeg = {
      if (!requireNamespace("jpeg", quietly = TRUE)) {
        cli::cli_abort(
          "Reading JPEG masks requires the {.pkg jpeg} package."
        )
      }
      jpeg::readJPEG(path)
    },
    cli::cli_abort(
      "Unsupported mask extension {.val {ext}} for {.path {path}}."
    )
  )
  if (length(dim(img)) == 3L) {
    img <- 0.2126 * img[, , 1] +
           0.7152 * img[, , 2] +
           0.0722 * img[, , 3]
  }
  img_dims <- as.integer(dim(img))
  if (!is.null(expected_dims)) {
    expected_dims <- as.integer(expected_dims)
    if (!identical(img_dims, expected_dims)) {
      cli::cli_abort(c(
        "Mask dimensions do not match {.arg expected_dims}.",
        "*" = "Loaded mask:    {img_dims[1]} x {img_dims[2]}",
        "*" = "Expected:       \\
               {expected_dims[1]} x {expected_dims[2]}"
      ))
    }
  }
  out <- as.vector(t(img) > threshold)
  if (isTRUE(invert)) out <- !out
  out
}

#' Internal: ellipse-mask helper
#'
#' Builds a logical matrix of pixels inside an ellipse with the
#' given centre and half-axes (all in 0-1 image fractions).
#' Centre = (row, col) in the 0-1 range; half-axes are in 0-1 of
#' the corresponding image side.
#'
#' @keywords internal
#' @noRd
ellipse_mask <- function(img_dims, c_row, c_col,
                         half_height, half_width) {
  nr <- img_dims[1L]; nc <- img_dims[2L]
  rr <- (row(matrix(0, nr, nc)) - 1L) / (nr - 1L) - c_row
  cc <- (col(matrix(0, nr, nc)) - 1L) / (nc - 1L) - c_col
  (rr / half_height)^2 + (cc / half_width)^2 <= 1
}
