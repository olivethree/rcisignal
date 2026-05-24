#' Per-pixel agreement heatmap for a producer signal matrix
#'
#' @description
#' Visualises where producers in a single condition agree on the
#' direction of signal. For each pixel, computes a one-sample
#' t-statistic against zero across producers
#' (`mean / (sd / sqrt(N))`), then displays the resulting map with
#' a diverging color palette (positive = agreement on positive
#' signal, negative = agreement on negative signal, zero = no
#' agreement). Saturation of the color is the magnitude of the
#' agreement, *not* the value of the group-mean CI.
#'
#' Use this to answer "where do producers consistently *agree* the
#' target trait lives in the face?". Pair with the group-mean CI
#' image (raw mask or rendered) to see *direction* and *agreement*
#' side by side.
#'
#' @details
#' This is structurally a one-sample t-map (vs zero); pixels where
#' producers' contributions are large *and consistent in sign* get
#' high |t|, pixels where contributions are random get t near zero.
#' Cluster-permutation inference would normally accompany this for
#' formal pixel-level FWER control between conditions
#' ([rel_cluster_test()]); the agreement map is the **descriptive**
#' counterpart for a single condition.
#'
#' @section Reading the plot:
#' Two palettes are available; pick by what question you are
#' asking the data.
#'
#' **`palette = "diverging"` (default).** Encodes sign and
#' magnitude together. Both deep red and deep blue indicate
#' **strong** agreement among producers; only the **direction**
#' differs. "No agreement" is the neutral color (white), not red.
#' * **Hue** encodes the sign of the per-pixel one-sample `t`.
#'   Blue = producers consistently *add* to the base at that pixel
#'   (positive agreement, producers chose noise that lightens the
#'   region); red = consistently *subtract* (negative agreement,
#'   producers chose noise that darkens the region).
#' * **Saturation** encodes `|t|`. Deep color at either end means
#'   strong, consistent agreement; pale color means weak or
#'   inconsistent. The colorbar on the right reads in `t` units.
#' * **`zlim`** is symmetric around zero by default so the neutral
#'   color aligns with `t = 0`. Pass `zlim = c(-z, z)` to fix the
#'   scale across panels for direct comparison.
#'
#' **`palette = "fire"`.** Encodes `|t|` only on a single-hue ramp
#' (pale yellow at zero -> deep red at large `|t|`). Use when the
#' question is *where* producers have a consistent opinion and the
#' direction is not needed. The `"fire"` view **discards sign by
#' design**; it cannot distinguish "producers consistently added"
#' from "producers consistently subtracted". To recover direction at
#' any region of interest, view the same data with
#' `palette = "diverging"` or pair with [plot_ci_overlay()] of the
#' group-mean CI.
#' * **Hue intensity** encodes `|t|`. Pale yellow / near-white at
#'   low `|t|` (so the underlying base face shows through low-
#'   agreement regions); orange at moderate `|t|`; deep red at
#'   large `|t|`. The colorbar reads in `|t|` units.
#' * **`zlim`** defaults to `c(0, max(|t|))` and is asymmetric.
#'
#' **Common to both palettes.**
#' * **`threshold`** clips color to the neutral end below
#'   `|t| < threshold`, making strong-agreement clusters stand out.
#'   This is descriptive only; it does not provide FWER control. For
#'   inferential pixel significance, use [agreement_map_test()] and
#'   render its result directly via
#'   `plot(agreement_map_test(...))`, or overlay the contours via
#'   [plot_ci_overlay()].
#' * **`base_image`** composites the heatmap on top of a grayscale
#'   base face so anatomical context shows through. Out-of-mask and
#'   subthreshold pixels render fully transparent; the per-pixel
#'   opacity scales `|t| / zlim_max` up to `alpha_max`. Works for
#'   both palettes. The color bar still shows the full scale so
#'   magnitudes are readable off the rendered overlay.
#' * The diverging color convention (blue = positive,
#'   red = negative) matches [plot_ci_overlay()] and the
#'   cluster-test plots so the same group CI reads consistently
#'   across the package. The `"fire"` option is unique to this
#'   function; the CI-overlay and cluster-test plots need to show
#'   direction, so they do not provide a magnitude-only view.
#'
#' @param signal_matrix Pixels x participants raw mask (as returned
#'   by `ci_from_responses_*()` or `read_cis()` + `extract_signal()`).
#' @param img_dims Integer `c(nrow, ncol)`. If `NULL`, inferred from
#'   `attr(signal_matrix, "img_dims")` or from `sqrt(n_pixels)` if
#'   the latter is a whole number.
#' @param mask Optional logical vector of length
#'   `nrow(signal_matrix)` (column-major) restricting display to a
#'   region (e.g., `make_face_mask(img_dims, region = "eyes")`).
#'   Also accepts the output of [read_face_mask()] for PNG/JPEG
#'   masks. Pixels outside the mask render as `NA` (transparent).
#' @param threshold Optional positive numeric. When supplied,
#'   pixels with `|t| < threshold` are rendered in the neutral
#'   (white) color, making clusters of agreement stand out.
#'   Default `NULL` (full continuous map).
#' @param zlim Numeric `c(low, high)` for the color scale. For
#'   `palette = "diverging"` (default), defaults to
#'   `c(-max(|t|), max(|t|))` so the neutral color aligns with
#'   `t = 0`. For `palette = "fire"`, defaults to `c(0, max(|t|))`
#'   so pale yellow aligns with `|t| = 0`.
#' @param palette Character. `"diverging"` (default; positive =
#'   blue, negative = red, neutral = white) encodes sign in hue and
#'   magnitude in saturation. `"fire"` encodes `|t|` only on a
#'   single-hue ramp (pale yellow at zero -> deep red at large
#'   `|t|`); use this when the question is "where do producers have
#'   a consistent opinion" and direction is not needed. The `"fire"`
#'   view discards sign; pair with `palette = "diverging"` or with
#'   [plot_ci_overlay()] to recover direction at a region of
#'   interest.
#' @param base_image Optional. Either a numeric matrix
#'   (`nrow x ncol`, grayscale, values in 0-1) or a path to a
#'   PNG/JPEG file. When supplied, the t-map is composited on top of
#'   the grayscale base; out-of-mask and subthreshold pixels render
#'   fully transparent. When `NULL` (default), the map is drawn on a
#'   flat panel via `graphics::image()` (the historical behavior).
#' @param alpha_max Numeric in `[0, 1]`. Maximum opacity of the
#'   heatmap at the color-scale top (`zlim_max`) when `base_image`
#'   is supplied. Ignored otherwise. Default 0.7.
#' @param main Title.
#' @param ... Passed to `graphics::image()`.
#' @return Invisibly, a list with `t_map` (numeric vector of t values
#'   per pixel; always signed regardless of palette), `n` (producer
#'   count), `img_dims`, `mask` (if supplied), `zlim` (the color
#'   scale used), and `palette` (the palette name).
#' @seealso [plot_ci_overlay()] for the producer-mean counterpart
#'   (signed CI, optionally with FWE contours); [agreement_map_test()]
#'   for FWE-controlled significance, and its `plot()` method for a
#'   one-call agreement map with contours; [rel_cluster_test()] for
#'   inferential between-condition tests; [make_face_mask()] /
#'   [read_face_mask()] for the optional `mask`.
#' @export
#' @examples
#' \dontrun{
#' # Minimal call-signature demo with a synthetic input. The agreement
#' # map will look flat because the input is pure noise.
#' n_side <- 32L
#' n_pix  <- n_side * n_side
#' set.seed(1)
#' signal_matrix <- matrix(rnorm(n_pix * 20L), n_pix, 20L)
#' plot_agreement_map(signal_matrix, img_dims = c(n_side, n_side))
#' }
#'
#' \dontrun{
#' # Same function, richer input: simulate Brief-RC responses with a
#' # signal planted in the eye region, then look at the agreement map.
#' # Producers should consistently agree on the planted region.
#' sim <- simulate_briefrc_data(
#'   n_per_condition = 20, n_trials = 60, conditions = "target",
#'   signal_region = "eyes", signal_strength = "strong", seed = 1
#' )
#' cis <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
#' plot_agreement_map(cis$signal_matrix)
#' }
#'
#' \dontrun{
#' # Composite the agreement map on the base face for a single
#' # publication-grade figure. Works for both palettes; the
#' # "diverging" branch matches plot_ci_overlay()'s color mapping.
#' sim <- simulate_briefrc_data(
#'   n_per_condition = 20, n_trials = 60, conditions = "target",
#'   signal_region = "eyes", signal_strength = "strong", seed = 1
#' )
#' cis <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
#' plot_agreement_map(cis$signal_matrix,
#'                    base_image = sim$base_face,
#'                    threshold  = 2.0,
#'                    main       = "Agreement t-map over base face")
#' }
plot_agreement_map <- function(signal_matrix,
                               img_dims  = NULL,
                               mask      = NULL,
                               threshold = NULL,
                               zlim      = NULL,
                               palette   = c("diverging", "fire"),
                               base_image = NULL,
                               alpha_max  = 0.7,
                               main      = "Per-pixel producer agreement (t-map)",
                               ...) {
  if (!is.matrix(signal_matrix) || !is.numeric(signal_matrix)) {
    cli::cli_abort("{.arg signal_matrix} must be a numeric matrix.")
  }
  palette <- match.arg(palette)

  n_pix <- nrow(signal_matrix)
  if (is.null(img_dims)) {
    a <- attr(signal_matrix, "img_dims")
    if (!is.null(a)) {
      img_dims <- a
    } else {
      side <- sqrt(n_pix)
      if (side != as.integer(side)) {
        cli::cli_abort(
          "Cannot infer {.arg img_dims}; pass it explicitly."
        )
      }
      img_dims <- c(as.integer(side), as.integer(side))
    }
  }

  n   <- ncol(signal_matrix)
  m   <- rowMeans(signal_matrix)
  v   <- rowSums((signal_matrix - m)^2) / (n - 1L)
  se  <- sqrt(v / n)
  t_map <- m / se
  t_map[!is.finite(t_map)] <- 0

  res <- render_agreement_t_map(
    t_map      = t_map,
    img_dims   = img_dims,
    mask       = mask,
    threshold  = threshold,
    zlim       = zlim,
    palette    = palette,
    base_image = base_image,
    alpha_max  = alpha_max,
    main       = main,
    sub_n      = n,
    ...
  )
  invisible(list(t_map = t_map, n = n,
                 img_dims = img_dims, mask = mask,
                 zlim = res$zlim, palette = palette))
}

#' Render an agreement-style t-map to the active device
#'
#' Internal renderer shared by `plot_agreement_map()` and
#' `plot.rcisignal_rel_agreement_map_test()`. Takes a pre-computed
#' per-pixel t-vector rather than a producer x pixel matrix so the
#' inferential plot method can pass `observed_t` straight through
#' without re-running the one-sample t-statistics.
#'
#' @keywords internal
#' @noRd
render_agreement_t_map <- function(t_map,
                                   img_dims,
                                   mask       = NULL,
                                   threshold  = NULL,
                                   zlim       = NULL,
                                   palette    = c("diverging", "fire"),
                                   base_image = NULL,
                                   alpha_max  = 0.7,
                                   main       = NULL,
                                   sub_n      = NULL,
                                   ...) {
  palette <- match.arg(palette)
  img_dims <- as.integer(img_dims)
  n_pix <- length(t_map)

  display <- if (palette == "fire") abs(t_map) else t_map
  if (!is.null(mask)) {
    if (!is.logical(mask) || length(mask) != n_pix) {
      cli::cli_abort(
        "{.arg mask} must be a logical vector of length {n_pix}."
      )
    }
    display[!mask] <- NA
  }
  if (!is.null(threshold)) {
    display[abs(display) < threshold] <- 0
  }

  if (is.null(zlim)) {
    rng <- max(abs(display), na.rm = TRUE)
    if (!is.finite(rng) || rng == 0) rng <- 1
    zlim <- if (palette == "fire") c(0, rng) else c(-rng, rng)
  }

  col_vec <- if (palette == "diverging") {
    grDevices::hcl.colors(256L, "RdBu")
  } else {
    grDevices::hcl.colors(256L, "YlOrRd", rev = TRUE)
  }
  bar_label <- if (palette == "fire") {
    "|t-value| (one-sample vs 0)"
  } else {
    "t-value (one-sample vs 0)"
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mar = c(1, 1, 3, 6) + 0.1)
  mat <- matrix(display, img_dims[1L], img_dims[2L])

  if (is.null(base_image)) {
    graphics::image(
      seq_len(img_dims[2L]), seq_len(img_dims[1L]),
      t(mat[nrow(mat):1L, ]),
      col       = col_vec,
      zlim      = zlim,
      main      = main,
      axes      = FALSE,
      xlab      = "", ylab = "",
      asp       = img_dims[1L] / img_dims[2L],
      useRaster = TRUE,
      ...
    )
    graphics::box(col = "grey80", lwd = 0.5)
  } else {
    base_mat <- resolve_base_for_overlay(base_image)
    if (!identical(as.integer(dim(base_mat)), as.integer(img_dims))) {
      cli::cli_abort(c(
        "{.arg base_image} dimensions do not match {.arg img_dims}.",
        "*" = "base_image: {dim(base_mat)[1L]} x {dim(base_mat)[2L]}",
        "*" = "img_dims:   {img_dims[1L]} x {img_dims[2L]}"
      ))
    }
    zlim_max <- max(abs(zlim))
    if (!is.finite(zlim_max) || zlim_max == 0) zlim_max <- 1

    display_alpha <- display
    display_alpha[is.na(display_alpha)] <- 0
    alpha_raster <- pmin(abs(display_alpha) / zlim_max, 1) * alpha_max
    alpha_raster <- matrix(alpha_raster, img_dims[1L], img_dims[2L])

    if (palette == "fire") {
      mag_mat <- matrix(abs(display_alpha), img_dims[1L], img_dims[2L])
      fg_rgb  <- fire_rgb_array(mag_mat, max_mag = zlim_max)
    } else {
      signed_mat <- matrix(display_alpha, img_dims[1L], img_dims[2L])
      fg_rgb     <- diverging_rgb_array(signed_mat)
    }
    composed <- composite_rgb_over_gray(base_mat, fg_rgb, alpha_raster)

    graphics::plot.new()
    graphics::plot.window(xlim = c(0, img_dims[2L]),
                          ylim = c(0, img_dims[1L]),
                          asp = 1, xaxs = "i", yaxs = "i")
    graphics::rasterImage(composed, 0, 0, img_dims[2L], img_dims[1L],
                          interpolate = FALSE)
    graphics::box(col = "grey80", lwd = 0.5)
    if (!is.null(main)) {
      graphics::title(main = main, line = 1, cex.main = 1.0,
                      font.main = 1)
    }
  }

  add_color_bar(zlim, col_vec, label = bar_label)

  if (!is.null(sub_n)) {
    graphics::mtext(
      sprintf("N = %d producers,  %d x %d pixels%s",
              sub_n, img_dims[1L], img_dims[2L],
              if (!is.null(threshold))
                sprintf(",  thresholded at |t| > %.2f", threshold)
              else ""),
      side = 3, line = 0.3, cex = 0.85, col = "grey30"
    )
  }
  invisible(list(zlim = zlim, palette = palette))
}

#' Add a vertical color bar in the right margin of the active plot
#'
#' @keywords internal
#' @noRd
add_color_bar <- function(zlim, col, label = NULL,
                           bar_width_frac = 0.05) {
  usr <- graphics::par("usr")
  pin <- graphics::par("pin")
  # Place bar in the right margin
  x_left  <- usr[2] + (usr[2] - usr[1]) * 0.03
  x_right <- usr[2] + (usr[2] - usr[1]) * 0.07
  y_bot   <- usr[3]
  y_top   <- usr[4]
  n_col   <- length(col)
  ys      <- seq(y_bot, y_top, length.out = n_col + 1L)
  graphics::par(xpd = TRUE)
  on.exit(graphics::par(xpd = FALSE), add = TRUE)
  for (i in seq_len(n_col)) {
    graphics::rect(x_left, ys[i], x_right, ys[i + 1L],
                   col = col[i], border = NA)
  }
  graphics::rect(x_left, y_bot, x_right, y_top,
                 border = "grey60", lwd = 0.5)
  ticks <- pretty(zlim, n = 5L)
  ticks <- ticks[ticks >= zlim[1] & ticks <= zlim[2]]
  tick_y <- y_bot + (ticks - zlim[1]) / diff(zlim) * (y_top - y_bot)
  graphics::segments(x_right, tick_y,
                     x_right + (usr[2] - usr[1]) * 0.01,
                     tick_y, col = "grey60", lwd = 0.5)
  graphics::text(x_right + (usr[2] - usr[1]) * 0.015, tick_y,
                 labels = format(ticks), pos = 4, cex = 0.7,
                 col = "grey20")
  if (!is.null(label)) {
    graphics::text(x_right + (usr[2] - usr[1]) * 0.06,
                   (y_bot + y_top) / 2,
                   labels = label, srt = -90, cex = 0.75,
                   col = "grey30")
  }
}
