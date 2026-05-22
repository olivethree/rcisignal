#' Pairwise Euclidean distance matrix across multiple group-mean CIs
#'
#' @description
#' Renders a publication-ready Euclidean distance matrix across
#' two or more classification images (CIs). A stage-2 consumer:
#' accepts a [group_ci()] result directly, a named numeric matrix
#' of `pixels x n_groups`, or the same named-list-of-CIs format as
#' [plot_ci_correlogram()] (per-producer matrices reduced to group
#' means internally).
#'
#' Use this when the question is "how *far apart* are these CIs in
#' pixel space?" rather than "how do they covary?". Euclidean
#' distance does not share the positive-baseline issue Pearson `r`
#' has on base-subtracted CIs (see [rel_dissimilarity()] for the
#' two-condition bootstrap version of the same magnitude metric).
#'
#' @section Reading the plot:
#' Each cell shows the Euclidean distance between two group-mean
#' CIs computed across the pixels included in `mask` (all pixels
#' by default). Larger values mean the two CIs sit farther apart
#' in pixel space. Saturation encodes magnitude on a single-hue
#' sequential ramp: pale (or near-white) at distance ~ 0, deep /
#' dark at large distances. Numbers in each cell show the exact
#' distance to `value_digits` decimal places.
#'
#' The diagonal is identically 0 (a CI's distance from itself).
#' With `triangle = "full"` it is rendered as the lightest cell.
#' With `triangle = "upper"` or `triangle = "lower"` the
#' diagonal cells are blank and carry the CI labels instead,
#' freeing the figure margins for paper figures.
#'
#' @section Raw vs normalised distance:
#' `method = "raw"` (default) returns the absolute Euclidean
#' distance, in whatever units the CIs carry. Useful for comparing
#' contrasts computed on the same image at the same resolution.
#'
#' `method = "normalised"` divides by `sqrt(n_pixels_used)`,
#' producing a per-pixel root-mean-square distance that is
#' comparable across masks and resolutions. Use this when the
#' question is "is contrast A more separated than contrast B"
#' and the two contrasts were computed under different masks (or
#' you want a publication number whose order of magnitude does
#' not balloon with image size).
#'
#' Both matrices are always returned (see "Reading the result"
#' below); the `method` argument only controls which is rendered
#' in the figure and used for cell labels.
#'
#' @section Reading the result:
#' The function returns invisibly a list with:
#' * `$distance_matrix`: the n x n distance matrix actually
#'   rendered (raw or normalised per `method`). Rows and columns
#'   are named with the CI names.
#' * `$distance_raw`: the raw distance matrix (always; useful for
#'   downstream MDS / `hclust()` regardless of which the figure
#'   showed).
#' * `$method`, `$n_pixels_used`, `$mask`, `$palette`,
#'   `$triangle`, `$file`: bookkeeping.
#'
#' @param cis CIs to compare, in any of three forms:
#'   - a [group_ci()] result (an [rcisignal_group_ci] matrix);
#'   - a numeric matrix `n_pixels x n_groups` with named columns;
#'   - a named list of CIs (each element a vector of length
#'     `prod(img_dims)`, a single-column matrix, or a per-producer
#'     `signal_matrix` from `ci_from_responses_*()` which is
#'     reduced to a group mean internally).
#'   Names become the row / column / diagonal labels in the
#'   figure. At least two CIs are required.
#' @param img_dims Integer `c(nrow, ncol)`. If `NULL`, inferred
#'   from `attr(cis[[1]], "img_dims")` or from `sqrt(n_pixels)`.
#' @param mask One of `"none"` (default), `"face"`, `"upper_face"`,
#'   or `"lower_face"`. Restricts the pixels included in the
#'   distance computation via [make_face_mask()].
#' @param method One of `"raw"` (default; absolute Euclidean) or
#'   `"normalised"` (`raw / sqrt(n_pixels_used)`).
#' @param triangle One of `"full"` (default), `"upper"` (mask the
#'   lower triangle and diagonal; CI labels on the diagonal
#'   cells), or `"lower"`.
#' @param palette One of `"viridis"` (default), `"inferno"`,
#'   `"plasma"`, or `"rocket"`. All four are sequential,
#'   colorblind-friendly palettes appropriate for non-negative
#'   distances.
#' @param show_values Logical. If `TRUE` (default), render the
#'   distance values inside each visible cell.
#' @param value_digits Integer or `NULL`. Decimal places for the
#'   cell labels. When `NULL` (default): 2 for raw distance, 3
#'   for normalised distance.
#' @param main Optional plot title.
#' @param file Optional output path. If `NULL` (default), plots
#'   to the current open device. If a path ending in `.png` or
#'   `.pdf` (case-insensitive), saves at 600 dpi (PNG) or as
#'   vector PDF.
#' @param width,height Optional output dimensions in inches.
#'   Defaults size the canvas to be square and scale with the
#'   number of CIs.
#' @param ... Currently unused; reserved for future arguments.
#' @return Invisibly, a list with `distance_matrix`,
#'   `distance_raw`, `method`, `n_pixels_used`, `mask`, `palette`,
#'   `triangle`, and `file`. See "Reading the result".
#' @seealso [plot_ci_correlogram()] for the Pearson-r version of
#'   the same input format; [rel_dissimilarity()] for a
#'   two-condition Euclidean distance with bootstrap CI;
#'   [plot_ci_mds()] to project multiple CIs into a 2D MDS scatter
#'   using the same distance matrix.
#' @export
#' @examples
#' \dontrun{
#' # Minimal: synthetic CIs to see the function's call signature
#' # and inspect the output shape.
#' set.seed(1)
#' n_pix <- 32L * 32L
#' ci_list <- list(
#'   A = rnorm(n_pix),
#'   B = rnorm(n_pix) + 0.3,
#'   C = rnorm(n_pix) - 0.2,
#'   D = rnorm(n_pix) + 0.5
#' )
#' out <- plot_ci_distance_matrix(ci_list, img_dims = c(32L, 32L))
#'
#' # The distance matrix (named numeric matrix):
#' out$distance_matrix
#'
#' # Raw distance is always returned, even if `method = "normalised"`:
#' out$distance_raw
#' }
#'
#' \dontrun{
#' # Realistic: simulate three conditions with planted signals in
#' # different face regions, build CIs, then compare in distance space.
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
#' plot_ci_distance_matrix(
#'   list("Eyes" = cis_eyes$signal_matrix,
#'        "Mouth" = cis_mouth$signal_matrix,
#'        "Nose"  = cis_nose$signal_matrix),
#'   mask     = "face",
#'   method   = "normalised",   # comparable across masks / resolutions
#'   triangle = "upper",
#'   file     = "ci_distance_matrix.pdf"
#' )
#' }
plot_ci_distance_matrix <- function(cis,
                                    img_dims     = NULL,
                                    mask         = c("none", "face",
                                                     "upper_face",
                                                     "lower_face"),
                                    method       = c("raw",
                                                     "normalised"),
                                    triangle     = c("full", "upper",
                                                     "lower"),
                                    palette      = c("viridis", "inferno",
                                                     "plasma", "rocket"),
                                    show_values  = TRUE,
                                    value_digits = NULL,
                                    main         = NULL,
                                    file         = NULL,
                                    width        = NULL,
                                    height       = NULL,
                                    ...) {
  mask     <- match.arg(mask)
  method   <- match.arg(method)
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

  d_raw <- as.matrix(stats::dist(t(M_use), method = "euclidean"))
  rownames(d_raw) <- nm
  colnames(d_raw) <- nm

  d_show <- if (identical(method, "raw")) {
    d_raw
  } else {
    d_raw / sqrt(n_pix_used)
  }

  if (is.null(value_digits)) {
    value_digits <- if (identical(method, "raw")) 2L else 3L
  }

  pal_name <- switch(palette,
                     viridis = "Viridis",
                     inferno = "Inferno",
                     plasma  = "Plasma",
                     rocket  = "Rocket")
  col_vec <- grDevices::hcl.colors(256L, pal_name, rev = TRUE)

  display <- d_show
  if (identical(triangle, "upper")) {
    display[lower.tri(display, diag = TRUE)] <- NA_real_
  } else if (identical(triangle, "lower")) {
    display[upper.tri(display, diag = TRUE)] <- NA_real_
  }

  zlim_max <- max(d_show, na.rm = TRUE)
  if (!is.finite(zlim_max) || zlim_max == 0) zlim_max <- 1
  zlim <- c(0, zlim_max)

  n_cis <- nrow(d_raw)
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
    zlim      = zlim,
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
    sat_threshold <- 0.55 * zlim_max
    for (i in seq_len(n_cis)) {
      for (j in seq_len(n_cis)) {
        v <- display[i, j]
        if (!is.na(v) && i != j) {
          graphics::text(
            x = j, y = n_cis - i + 1L,
            labels = sprintf("%.*f", value_digits, v),
            col = if (v > sat_threshold) "white" else "grey20",
            cex = 0.85
          )
        }
      }
    }
  }

  if (!is.null(main)) {
    graphics::mtext(main, side = 3, line = top_mar - 1.2,
                    cex = 1.0, col = "grey15")
  }

  bar_label <- if (identical(method, "raw")) {
    "Euclidean distance"
  } else {
    "Euclidean / sqrt(n_pixels)"
  }
  add_colour_bar(zlim, col_vec, label = bar_label)

  invisible(list(
    distance_matrix = d_show,
    distance_raw    = d_raw,
    method          = method,
    n_pixels_used   = n_pix_used,
    mask            = mask,
    palette         = palette,
    triangle        = triangle,
    file            = file
  ))
}
