#' Pairwise Pearson-r correlogram across multiple group-mean CIs
#'
#' @description
#' Renders a publication-ready Pearson correlation matrix across
#' two or more classification images (CIs), with optional face
#' region masking and optional direct save to PNG (600 dpi) or
#' PDF.
#'
#' Designed for the common publication-figure task: given a small
#' set of group-mean CIs from different conditions, show how
#' similar they are pixel-wise. The function accepts CIs in three
#' shapes (column-major vectors, single-column matrices, or
#' per-producer `signal_matrix` returned by
#' `ci_from_responses_*()`); per-producer matrices are reduced to
#' group means via `rowMeans()` automatically.
#'
#' Diverging-palette convention matches the rest of the package:
#' positive `r` = blue, negative `r` = red, neutral colour at
#' zero. The colour scale is fixed at `c(-1, 1)` so panels are
#' directly comparable across runs and across paper figures.
#'
#' @section Reading the plot:
#' Each cell shows the Pearson correlation `r` between two
#' group-mean CIs computed across the pixels included in `mask`
#' (all pixels by default). Saturation encodes `|r|`; hue encodes
#' sign (blue = positive, red = negative). Numbers in each cell
#' show the exact `r` to `value_digits` decimal places.
#'
#' The diagonal is identically 1 by definition. With
#' `triangle = "full"` the diagonal is rendered as a saturated
#' blue cell with `+1.00` inside. With `triangle = "upper"` or
#' `triangle = "lower"` the diagonal cells are blank and carry
#' the CI labels instead, freeing the figure margins.
#'
#' @section Pearson r between CIs is not a clean similarity score:
#' Two base-subtracted CIs share image-domain spatial structure
#' (face shape, oval signal support) that pushes their
#' correlation above zero even when the underlying mental
#' representations are unrelated. A high absolute `r` does not
#' by itself license a similarity claim. Use this correlogram to
#' visualise the **relative** ordering of pair-wise correlations
#' (which pairs covary more than others), and benchmark any
#' absolute claim against a permutation null. For a
#' baseline-free magnitude summary of how different two
#' conditions are, see [rel_dissimilarity()] (Euclidean
#' distance with bootstrap CI). For pixel-level localisation
#' of *where* two conditions differ, see [rel_cluster_test()].
#'
#' @param cis A named list of CIs to correlate. Each element may
#'   be a numeric vector of length `prod(img_dims)` (a single
#'   group-mean CI in column-major order), a numeric matrix
#'   `n_pixels x 1`, or a numeric matrix `n_pixels x n_producers`
#'   (a per-producer `signal_matrix` from `ci_from_responses_*()`;
#'   the function takes `rowMeans()` to obtain the group mean).
#'   Names become the row / column / diagonal labels in the
#'   figure. At least two elements are required.
#' @param img_dims Integer `c(nrow, ncol)`. If `NULL`, inferred
#'   from `attr(cis[[1]], "img_dims")` (set by `ci_from_responses_*()`)
#'   or from `sqrt(n_pixels)` if that is a whole number.
#' @param mask One of `"none"` (default; correlate over all
#'   pixels), `"face"` (full-face oval from [make_face_mask()]),
#'   `"upper_face"` (top half of the oval), or `"lower_face"`
#'   (bottom half). Restricting to the face oval removes
#'   off-face noise from the correlation; sub-regions answer
#'   "do these CIs covary in the upper / lower face only".
#' @param triangle One of `"full"` (default; render the whole
#'   matrix), `"upper"` (mask the lower triangle and diagonal;
#'   put CI labels on the diagonal cells), or `"lower"` (mirror).
#' @param palette One of `"diverging"` (default; RdBu, blue =
#'   positive), `"diverging_puor"` (PuOr, purple = positive), or
#'   `"diverging_brbg"` (BrBG, green = positive). All three are
#'   colorblind-friendly diverging palettes appropriate for
#'   correlation matrices.
#' @param show_values Logical. If `TRUE` (default), render the
#'   correlation values inside each visible cell.
#' @param value_digits Integer. Decimal places for the cell-value
#'   labels. Default `2L` (e.g., `"+0.67"`).
#' @param main Optional plot title.
#' @param file Optional output path. If `NULL` (default), plots
#'   to the current open device. If a path ending in `.png` or
#'   `.pdf` (case-insensitive), opens the corresponding device
#'   and writes the figure to disk: PNG at 600 dpi (raster), PDF
#'   as vector. Anything else aborts.
#' @param width,height Optional output dimensions in inches.
#'   Default sizes the canvas to be square and scale with the
#'   number of CIs: `max(5.5, 0.55 * n_cis + 4.5)` inches per
#'   side. Override either to control the final figure size.
#' @param ... Currently unused; reserved for future arguments.
#' @return Invisibly, a list with `correlation_matrix` (the n x n
#'   named Pearson-r matrix actually plotted), `n_pixels_used`,
#'   `mask`, `palette`, `triangle`, and `file` (the path written
#'   to, or `NULL`).
#' @seealso [rel_dissimilarity()] for a Euclidean-distance
#'   (baseline-free) magnitude summary between two conditions;
#'   [rel_cluster_test()] for pixel-level inference on where
#'   two conditions differ; [make_face_mask()] for the masks
#'   used by `mask = `.
#' @export
#' @examples
#' \dontrun{
#' # Minimal call-signature demo with synthetic vectors.
#' n_side <- 32L
#' n_pix  <- n_side * n_side
#' set.seed(1)
#' ci_list <- list(
#'   A = rnorm(n_pix),
#'   B = rnorm(n_pix),
#'   C = rnorm(n_pix)
#' )
#' plot_ci_correlogram(ci_list, img_dims = c(n_side, n_side))
#' }
#'
#' \dontrun{
#' # Realistic input: simulate three conditions with planted signals
#' # in different face regions, then compare their CIs.
#' sim_eyes  <- simulate_briefrc_data(
#'   n_per_condition = 20, n_trials = 60, conditions = "eyes",
#'   signal_region = "eyes", signal_strength = "strong", seed = 1
#' )
#' sim_mouth <- simulate_briefrc_data(
#'   n_per_condition = 20, n_trials = 60, conditions = "mouth",
#'   signal_region = "mouth", signal_strength = "strong", seed = 2
#' )
#' sim_nose  <- simulate_briefrc_data(
#'   n_per_condition = 20, n_trials = 60, conditions = "nose",
#'   signal_region = "nose", signal_strength = "strong", seed = 3
#' )
#' cis_eyes  <- ci_from_responses_briefrc(sim_eyes$data,
#'                                        noise_matrix = sim_eyes$noise_matrix)
#' cis_mouth <- ci_from_responses_briefrc(sim_mouth$data,
#'                                        noise_matrix = sim_mouth$noise_matrix)
#' cis_nose  <- ci_from_responses_briefrc(sim_nose$data,
#'                                        noise_matrix = sim_nose$noise_matrix)
#'
#' plot_ci_correlogram(
#'   list("Eyes" = cis_eyes$signal_matrix,
#'        "Mouth" = cis_mouth$signal_matrix,
#'        "Nose"  = cis_nose$signal_matrix),
#'   mask     = "face",
#'   triangle = "upper",
#'   file     = "ci_correlogram.pdf"
#' )
#' }
plot_ci_correlogram <- function(cis,
                                img_dims     = NULL,
                                mask         = c("none", "face",
                                                 "upper_face",
                                                 "lower_face"),
                                triangle     = c("full", "upper",
                                                 "lower"),
                                palette      = c("diverging",
                                                 "diverging_puor",
                                                 "diverging_brbg"),
                                show_values  = TRUE,
                                value_digits = 2L,
                                main         = NULL,
                                file         = NULL,
                                width        = NULL,
                                height       = NULL,
                                ...) {
  mask     <- match.arg(mask)
  triangle <- match.arg(triangle)
  palette  <- match.arg(palette)

  prep      <- prepare_ci_matrix(cis, img_dims, min_cis = 2L)
  M         <- prep$M
  img_dims  <- prep$img_dims
  n_pix     <- prep$n_pix
  nm        <- prep$names

  mask_vec   <- resolve_ci_mask(mask, img_dims, n_pix)
  M_use      <- M[mask_vec, , drop = FALSE]
  n_pix_used <- sum(mask_vec)

  cor_mat <- stats::cor(M_use, method = "pearson")
  if (any(!is.finite(cor_mat))) {
    cli::cli_abort(c(
      "Correlation matrix contains non-finite values.",
      "i" = "Check that no CI has zero variance over the included pixels."
    ))
  }

  pal_name <- switch(palette,
                     diverging       = "RdBu",
                     diverging_puor  = "PuOr",
                     diverging_brbg  = "BrBG")
  col_vec <- grDevices::hcl.colors(256L, pal_name)

  display <- cor_mat
  if (identical(triangle, "upper")) {
    display[lower.tri(display, diag = TRUE)] <- NA_real_
  } else if (identical(triangle, "lower")) {
    display[upper.tri(display, diag = TRUE)] <- NA_real_
  }

  n_cis <- nrow(cor_mat)
  default_dim_in <- max(5.5, 0.55 * n_cis + 4.5)
  W <- if (is.null(width))  default_dim_in else width
  H <- if (is.null(height)) default_dim_in else height

  device_opened <- FALSE
  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (identical(ext, "png")) {
      grDevices::png(file, width = W, height = H,
                     units = "in", res = 600)
      device_opened <- TRUE
    } else if (identical(ext, "pdf")) {
      grDevices::pdf(file, width = W, height = H)
      device_opened <- TRUE
    } else {
      cli::cli_abort(c(
        "{.arg file} must end in {.val .png} or {.val .pdf}.",
        "i" = "Got extension {.val {ext}}."
      ))
    }
  }
  if (device_opened) on.exit(grDevices::dev.off(), add = TRUE)

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  top_mar <- if (identical(triangle, "full")) {
    max(5, max(nchar(nm)) * 0.55)
  } else {
    1.5
  }
  left_mar <- if (identical(triangle, "full")) {
    max(5, max(nchar(nm)) * 0.55)
  } else {
    1
  }
  graphics::par(mar = c(0.5, left_mar, top_mar, 7) + 0.1)

  graphics::image(
    seq_len(n_cis), seq_len(n_cis),
    t(display[n_cis:1L, , drop = FALSE]),
    col       = col_vec,
    zlim      = c(-1, 1),
    axes      = FALSE, xlab = "", ylab = "",
    asp       = 1, useRaster = TRUE
  )

  if (identical(triangle, "full")) {
    graphics::axis(3, at = seq_len(n_cis), labels = nm,
                   las = 2, tick = FALSE, line = -0.5,
                   cex.axis = 0.85)
    graphics::axis(2, at = rev(seq_len(n_cis)), labels = nm,
                   las = 1, tick = FALSE, line = -0.5,
                   cex.axis = 0.85)
  } else {
    for (i in seq_len(n_cis)) {
      graphics::text(
        x = i, y = n_cis - i + 1L,
        labels = nm[i], cex = 0.9, col = "grey15", font = 2
      )
    }
  }

  if (isTRUE(show_values)) {
    for (i in seq_len(n_cis)) {
      for (j in seq_len(n_cis)) {
        v <- display[i, j]
        if (!is.na(v) && i != j) {
          graphics::text(
            x = j, y = n_cis - i + 1L,
            labels = sprintf("%+.*f", value_digits, v),
            col = if (abs(v) > 0.6) "white" else "grey20",
            cex = 0.85
          )
        } else if (!is.na(v) && i == j && identical(triangle, "full")) {
          graphics::text(
            x = j, y = n_cis - i + 1L,
            labels = sprintf("%+.*f", value_digits, v),
            col = "white", cex = 0.85
          )
        }
      }
    }
  }

  if (!is.null(main)) {
    graphics::mtext(main, side = 3, line = top_mar - 1.2,
                    cex = 1.0, col = "grey15")
  }

  add_colour_bar(c(-1, 1), col_vec, label = "Pearson r")

  invisible(list(
    correlation_matrix = cor_mat,
    n_pixels_used      = n_pix_used,
    mask               = mask,
    palette            = palette,
    triangle           = triangle,
    file               = file
  ))
}
