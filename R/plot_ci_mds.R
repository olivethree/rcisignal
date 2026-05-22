#' Multidimensional-scaling (MDS) projection of multiple CIs
#'
#' @description
#' Projects two or more classification images (CIs) into a low-
#' dimensional Euclidean scatter using classical MDS
#' (`stats::cmdscale`). Distances between points in the scatter
#' reproduce the original Euclidean distances between CIs as
#' faithfully as the chosen number of dimensions allows.
#'
#' Use this to see at a glance which condition CIs cluster
#' together and which sit far apart, especially when you have
#' more than a handful of conditions to compare side by side.
#'
#' @section Auto-selecting the number of dimensions:
#' By default (`k = "auto"`), the function fits classical MDS at
#' every dimensionality from 2 up to `k_max` and picks the
#' **smallest** `k` whose Kruskal stress-1 against the original
#' Euclidean distances reaches the "good" threshold (default
#' `stress_threshold = 0.05`).
#'
#' Kruskal's (1964) interpretive bands for stress-1:
#'
#' * `0.025` excellent (the 2D / kD map is essentially exact)
#' * `0.05`  good      (small distortions; safe to interpret)
#' * `0.10`  fair      (interpret carefully; check the trace)
#' * `0.20`  poor      (the projection is hiding more than it shows)
#' * `> 0.20` very poor
#'
#' When the selected `k > 2`, the figure becomes a grid of all
#' `choose(k, 2)` pairwise dimension panels (Dim 1 vs 2, Dim 1
#' vs 3, Dim 2 vs 3, ...) so no information is hidden by a
#' premature flattening to 2D.
#'
#' The companion goodness-of-fit metric (Mardia's GOF1,
#' cumulative eigenvalue variance) is reported alongside stress
#' in `$variance_pct_by_k`. Values near 100% mean the kD map
#' captures essentially all of the positive eigenmass; below
#' ~60% the kD map is hiding more than it shows. Stress is
#' usually the more interpretable of the two.
#'
#' @section Theory-driven dimensionality:
#' Real research often has a theoretical reason to fix `k` to a
#' specific number (typically 2, for a paper-style scatter that
#' matches a two-axis hypothesis like warmth/dominance, or a
#' single panel for visual comparison with prior work). Pass an
#' integer `k` to bypass auto-selection:
#'
#' ```r
#' plot_ci_mds(ci_list, mask = "face", k = 2L)
#' ```
#'
#' The function still computes stress at every `k` in
#' `[2, k_max]` and exposes the trace via `$stress_by_k`. When
#' the forced dimensionality has high stress, report
#' `out$stress_by_k` alongside the 2D figure in the paper so
#' readers see what the theory-driven projection is hiding;
#' interpret point positions in *relative* terms (which
#' conditions cluster together) rather than as absolute
#' distances.
#'
#' Force a higher specific `k` (e.g., `k = 3L`) when a richer
#' projection is theoretically motivated; the function will
#' render the corresponding `choose(k, 2)`-panel grid.
#'
#' @section Reading the plot:
#' Each panel is a scatter of CIs in two MDS dimensions. Points
#' close together represent CIs that are similar in pixel space;
#' distant points represent dissimilar CIs. With `groups`,
#' points colour by category; with `shapes`, points carry
#' distinct marker shapes by another category. The aspect ratio
#' is fixed (`asp = 1`) so visual distances faithfully reflect
#' the MDS distances.
#'
#' The line at the top of the figure reports the selected `k`,
#' the Kruskal stress at that `k` and its interpretive band,
#' and the per-axis variance percentages.
#'
#' @section Reading the result:
#' The returned object (class `rcisignal_mds`) contains:
#'
#' * `$mds_points`: **the coordinates of each CI in the Euclidean
#'   MDS space**. An `n_cis x k_selected` named numeric matrix
#'   (rows = CI names, columns = orthogonal MDS axes). Pull these
#'   directly to compute custom downstream analyses (e.g.
#'   distances between specific CIs in the reduced space, your
#'   own scatter, clustering).
#' * `$distance_matrix`: the pairwise distances the projection
#'   was fit to (raw or normalised per `distance`).
#' * `$stress_by_k`, `$variance_pct_by_k`: per-`k` Kruskal stress
#'   and cumulative Mardia GOF1. Inspect these to audit the
#'   auto-selected `k`.
#' * `$k_selected`, `$stress_1`, `$stress_band`,
#'   `$variance_pct`: summaries at the selected `k`.
#' * `$panel_pairs`: the `2 x choose(k, 2)` matrix of dimension
#'   pairs the figure rendered.
#'
#' Call `print()` or `summary()` on the returned object for a
#' one-screen human-readable view.
#'
#' @param cis CIs to project, in any of three forms:
#'   - a [group_ci()] result (an [rcisignal_group_ci] matrix);
#'   - a numeric matrix `n_pixels x n_groups` with named columns;
#'   - a named list of CIs (each element a vector of length
#'     `prod(img_dims)`, a single-column matrix, or a per-producer
#'     `signal_matrix` from `ci_from_responses_*()` which is
#'     reduced to a group mean internally).
#'   At least three CIs are required for a meaningful 2D
#'   projection.
#' @param img_dims Integer `c(nrow, ncol)`. If `NULL`, inferred
#'   from `attr(cis[[1]], "img_dims")` or from `sqrt(n_pixels)`.
#' @param mask One of `"none"` (default), `"face"`, `"upper_face"`,
#'   or `"lower_face"`. Restricts the pixels included in the
#'   distance computation via [make_face_mask()].
#' @param distance One of `"euclidean_raw"` (default; absolute
#'   Euclidean distance) or `"euclidean_normalised"`
#'   (`raw / sqrt(n_pixels_used)`). The MDS projection is fit to
#'   whichever distance is selected.
#' @param k Either `"auto"` (default; pick the smallest `k`
#'   whose Kruskal stress reaches `stress_threshold`) or an
#'   integer (forces that many dimensions; bypasses auto-
#'   selection). Use `k = 2L` for a single paper-figure panel.
#' @param stress_threshold Numeric. The Kruskal stress level the
#'   auto-selector treats as "good enough". Default `0.05`
#'   (Kruskal's "good" band). Lower it to `0.025` for "excellent"
#'   fidelity; raise to `0.10` for "fair" if the data are too
#'   noisy to reach 0.05 at small `k`.
#' @param k_max Integer. Maximum dimensionality the auto-selector
#'   will try. Default `4L`, capping the multi-panel grid at
#'   `choose(4, 2) = 6` pairs. Silently capped to
#'   `length(cis) - 1L` (the maximum embeddable dimension).
#' @param groups Optional named character vector mapping CI names
#'   to a categorical group label. Used to colour the points.
#'   Names must match `names(cis)` exactly.
#' @param shapes Optional named character vector with the same
#'   shape as `groups`, mapping CI names to a second categorical
#'   label used for point shapes (`pch`).
#' @param point_size,point_alpha Numeric. Point cex and alpha.
#' @param label Logical. Draw CI names above each point. Default
#'   `TRUE`.
#' @param label_cex Numeric. Cex for point labels.
#' @param show_axes_at_zero Logical. Draw faint zero-reference
#'   lines on each panel. Default `TRUE`.
#' @param show_gof Logical. Render the GOF summary above the
#'   figure (selected `k`, stress + band, per-axis variance).
#'   Default `TRUE`.
#' @param main Optional plot title.
#' @param file Optional output path. If `NULL` (default), plots
#'   to the current open device. If a path ending in `.png` or
#'   `.pdf` (case-insensitive), saves at 600 dpi (PNG) or as
#'   vector PDF.
#' @param width,height Optional output dimensions in inches.
#'   Defaults scale with the number of CIs and the number of
#'   panels.
#' @param ... Currently unused; reserved for future arguments.
#' @return Invisibly, an object of class `rcisignal_mds`. See
#'   "Reading the result".
#' @seealso [plot_ci_distance_matrix()] for the all-vs-all
#'   distance matrix that this function projects;
#'   [plot_ci_correlogram()] for the Pearson-r version of the
#'   same multi-CI summary; [rel_dissimilarity()] for a
#'   two-condition Euclidean distance with bootstrap CI.
#' @references
#' Kruskal, J. B. (1964). Multidimensional scaling by optimizing
#'   goodness of fit to a nonmetric hypothesis. *Psychometrika*,
#'   29(1), 1-27.
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
#'   D = rnorm(n_pix) + 0.5,
#'   E = rnorm(n_pix)
#' )
#'
#' # Auto-selects the smallest k reaching "good" stress (<= 0.05).
#' out <- plot_ci_mds(ci_list, img_dims = c(32L, 32L))
#'
#' # Human-readable view of the dimensionality selection trace:
#' print(out)
#'
#' # Extract the coordinates of each CI in the Euclidean MDS space.
#' out$mds_points
#'
#' # Distance between CI "A" and CI "C" in the reduced space:
#' sqrt(sum((out$mds_points["A", ] - out$mds_points["C", ])^2))
#'
#' # Per-k stress and variance traces (was the auto-choice sensible?):
#' out$stress_by_k
#' out$variance_pct_by_k
#' }
#'
#' \dontrun{
#' # Realistic: simulate four conditions with planted signals,
#' # build CIs, then compare in MDS space.
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
#' sim_flat  <- simulate_briefrc_data(
#'   n_per_condition = 20, n_trials = 60, conditions = "control",
#'   signal_region = NULL, seed = 4
#' )
#'
#' cis_eyes  <- ci_from_responses_briefrc(sim_eyes$data,
#'                                        noise_matrix = sim_eyes$noise_matrix)
#' cis_mouth <- ci_from_responses_briefrc(sim_mouth$data,
#'                                        noise_matrix = sim_mouth$noise_matrix)
#' cis_nose  <- ci_from_responses_briefrc(sim_nose$data,
#'                                        noise_matrix = sim_nose$noise_matrix)
#' cis_flat  <- ci_from_responses_briefrc(sim_flat$data,
#'                                        noise_matrix = sim_flat$noise_matrix)
#'
#' ci_list <- list(
#'   "Eyes"    = cis_eyes$signal_matrix,
#'   "Mouth"   = cis_mouth$signal_matrix,
#'   "Nose"    = cis_nose$signal_matrix,
#'   "Control" = cis_flat$signal_matrix
#' )
#'
#' out <- plot_ci_mds(ci_list, mask = "face")
#'
#' # Force a single 2D paper figure once you have audited fidelity:
#' plot_ci_mds(ci_list, mask = "face", k = 2L,
#'             file = "fig_mds.pdf")
#'
#' # Add a grouping variable for point colour:
#' plot_ci_mds(
#'   ci_list, mask = "face",
#'   groups = c(Eyes = "feature", Mouth = "feature",
#'              Nose = "feature", Control = "control")
#' )
#' }
plot_ci_mds <- function(cis,
                        img_dims          = NULL,
                        mask              = c("none", "face",
                                              "upper_face",
                                              "lower_face"),
                        distance          = c("euclidean_raw",
                                              "euclidean_normalised"),
                        k                 = "auto",
                        stress_threshold  = 0.05,
                        k_max             = 4L,
                        groups            = NULL,
                        shapes            = NULL,
                        point_size        = 2.5,
                        point_alpha       = 0.85,
                        label             = TRUE,
                        label_cex         = 0.85,
                        show_axes_at_zero = TRUE,
                        show_gof          = TRUE,
                        main              = NULL,
                        file              = NULL,
                        width             = NULL,
                        height            = NULL,
                        ...) {
  mask     <- match.arg(mask)
  distance <- match.arg(distance)

  prep      <- prepare_ci_matrix(cis, img_dims, min_cis = 3L)
  M         <- prep$M
  img_dims  <- prep$img_dims
  n_pix     <- prep$n_pix
  nm        <- prep$names
  n_cis     <- length(nm)

  if (!is.numeric(stress_threshold) ||
      stress_threshold <= 0 || stress_threshold > 1) {
    cli::cli_abort(
      "{.arg stress_threshold} must be a positive number in (0, 1]."
    )
  }
  if (!is.numeric(k_max) || length(k_max) != 1L || k_max < 2L) {
    cli::cli_abort("{.arg k_max} must be an integer >= 2.")
  }
  k_max <- as.integer(k_max)
  k_max <- min(k_max, n_cis - 1L)

  is_auto <- identical(k, "auto") || identical(k, "Auto") ||
             identical(k, "AUTO")
  if (!is_auto) {
    if (!is.numeric(k) || length(k) != 1L) {
      cli::cli_abort(
        "{.arg k} must be {.val \"auto\"} or a single integer."
      )
    }
    k_user <- as.integer(k)
    if (k_user < 2L || k_user > n_cis - 1L) {
      cli::cli_abort(c(
        "{.arg k} = {k_user} is out of range.",
        "i" = "Must be in [2, {n_cis - 1L}] for {n_cis} CIs."
      ))
    }
  }

  mask_vec   <- resolve_ci_mask(mask, img_dims, n_pix)
  M_use      <- M[mask_vec, , drop = FALSE]
  n_pix_used <- sum(mask_vec)

  d_raw <- as.matrix(stats::dist(t(M_use), method = "euclidean"))
  rownames(d_raw) <- nm
  colnames(d_raw) <- nm
  d_used <- if (identical(distance, "euclidean_raw")) {
    d_raw
  } else {
    d_raw / sqrt(n_pix_used)
  }

  k_grid <- 2L:k_max
  mds_fits <- lapply(k_grid, function(k_) {
    stats::cmdscale(stats::as.dist(d_used), k = k_, eig = TRUE)
  })

  drop_idx <- vapply(seq_along(mds_fits), function(i) {
    ncol(mds_fits[[i]]$points) < k_grid[i]
  }, logical(1L))
  if (all(drop_idx)) {
    cli::cli_abort(c(
      "Distance matrix is rank-deficient; classical MDS cannot \\
       produce a 2D embedding.",
      "i" = "Check for duplicate or zero-variance CIs in {.arg cis}."
    ))
  }
  k_grid   <- k_grid[!drop_idx]
  mds_fits <- mds_fits[!drop_idx]

  stress_by_k <- vapply(mds_fits, function(fit) {
    mds_d <- as.matrix(stats::dist(fit$points))
    upper <- upper.tri(d_used)
    sqrt(sum((d_used[upper] - mds_d[upper])^2) /
         sum(d_used[upper]^2))
  }, numeric(1L))
  names(stress_by_k) <- as.character(k_grid)

  abs_eig_full <- abs(mds_fits[[length(mds_fits)]]$eig)
  variance_pct_by_k <- vapply(k_grid, function(k_) {
    sum(abs_eig_full[seq_len(k_)]) / sum(abs_eig_full) * 100
  }, numeric(1L))
  names(variance_pct_by_k) <- as.character(k_grid)

  if (is_auto) {
    below <- which(stress_by_k <= stress_threshold)
    if (length(below)) {
      k_use <- as.integer(names(stress_by_k)[below[1L]])
    } else {
      k_use <- max(k_grid)
      best_k     <- as.integer(names(stress_by_k)[which.min(stress_by_k)])
      best_str   <- min(stress_by_k)
      cli::cli_warn(c(
        "No tested dimensionality reached stress <= {stress_threshold}.",
        "i" = paste0(
          "Best stress in [2, ", max(k_grid),
          "] = ", sprintf("%.3f", best_str),
          " at k = ", best_k, "."
        ),
        "i" = paste0(
          "Falling back to k = ", k_use,
          ". For a theory-driven 2D figure, pass {.code k = 2L} ",
          "explicitly and report {.code out$stress_by_k} alongside ",
          "the figure. Alternatively raise {.arg k_max}."
        )
      ))
    }
  } else {
    k_use <- k_user
    if (!(as.character(k_use) %in% names(stress_by_k))) {
      cli::cli_abort(
        "{.arg k} = {k_use} could not be embedded (rank-deficient)."
      )
    }
  }

  mds          <- mds_fits[[which(k_grid == k_use)]]
  mds_points   <- mds$points
  rownames(mds_points) <- nm
  colnames(mds_points) <- paste0("dim", seq_len(k_use))

  stress_1 <- unname(stress_by_k[as.character(k_use)])
  abs_eig_sel  <- abs(mds$eig)
  variance_pct <- abs_eig_sel[seq_len(k_use)] / sum(abs_eig_sel) * 100
  names(variance_pct) <- paste0("dim", seq_len(k_use))

  band <- if (stress_1 <= 0.025) "excellent"
          else if (stress_1 <= 0.05)  "good"
          else if (stress_1 <= 0.10)  "fair"
          else if (stress_1 <= 0.20)  "poor"
          else                        "very poor"

  panel_pairs <- if (k_use == 2L) {
    matrix(c(1L, 2L), nrow = 2L)
  } else {
    utils::combn(k_use, 2L)
  }
  n_panels <- ncol(panel_pairs)

  groups <- validate_categorical_arg(groups, nm, "groups")
  shapes <- validate_categorical_arg(shapes, nm, "shapes")

  if (!is.null(groups)) {
    grp_lev   <- sort(unique(groups))
    grp_cols  <- grDevices::hcl.colors(max(length(grp_lev), 2L),
                                       "Dark 3")[seq_along(grp_lev)]
    point_col <- grp_cols[match(groups, grp_lev)]
  } else {
    grp_lev   <- NULL
    grp_cols  <- "#404040"
    point_col <- rep("#404040", n_cis)
  }
  if (!is.null(shapes)) {
    shp_lev   <- sort(unique(shapes))
    shp_pchs  <- c(16, 17, 18, 15, 8, 7, 5, 4, 3)
    if (length(shp_lev) > length(shp_pchs)) {
      shp_pchs <- rep_len(shp_pchs, length(shp_lev))
    }
    point_pch <- shp_pchs[match(shapes, shp_lev)]
  } else {
    shp_lev   <- NULL
    shp_pchs  <- 16
    point_pch <- rep(16L, n_cis)
  }

  if (n_panels == 1L) {
    default_w <- max(6, 0.4 * n_cis + 4)
    default_h <- max(5.5, 0.35 * n_cis + 4)
  } else {
    nrow_grid <- ceiling(sqrt(n_panels))
    ncol_grid <- ceiling(n_panels / nrow_grid)
    default_w <- ncol_grid * 3.8 + 1.0
    default_h <- nrow_grid * 3.8 + 1.0
  }
  W <- if (is.null(width))  default_w else width
  H <- if (is.null(height)) default_h else height

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

  if (n_panels == 1L) {
    graphics::par(mar = c(4, 4.2, 3.4, 1) + 0.1)
  } else {
    nrow_grid <- ceiling(sqrt(n_panels))
    ncol_grid <- ceiling(n_panels / nrow_grid)
    graphics::par(mfrow = c(nrow_grid, ncol_grid),
                  oma   = c(0, 0, 3.2, 0),
                  mar   = c(3.6, 3.8, 1.6, 0.6) + 0.1)
  }

  point_col_alpha <- vapply(point_col, function(co) {
    grDevices::adjustcolor(co, alpha.f = point_alpha)
  }, character(1L))

  for (p in seq_len(n_panels)) {
    dx <- panel_pairs[1L, p]
    dy <- panel_pairs[2L, p]
    xs <- mds_points[, dx]
    ys <- mds_points[, dy]
    rng_x <- range(xs); rng_y <- range(ys)
    pad_x <- max(diff(rng_x) * 0.12, 1e-12)
    pad_y <- max(diff(rng_y) * 0.12, 1e-12)

    graphics::plot.new()
    graphics::plot.window(
      xlim = rng_x + c(-pad_x, pad_x),
      ylim = rng_y + c(-pad_y, pad_y),
      asp  = 1
    )
    if (isTRUE(show_axes_at_zero)) {
      graphics::abline(h = 0, v = 0, col = "grey85",
                       lwd = 0.5, lty = 1)
    }
    graphics::axis(1, col = "grey50", col.axis = "grey20",
                   cex.axis = 0.78, lwd = 0.7)
    graphics::axis(2, col = "grey50", col.axis = "grey20",
                   cex.axis = 0.78, lwd = 0.7, las = 1)
    graphics::box(col = "grey80", lwd = 0.5)
    graphics::title(xlab = sprintf("MDS dim %d", dx),
                    ylab = sprintf("MDS dim %d", dy),
                    line = 2.3, cex.lab = 0.85, col.lab = "grey20")

    graphics::points(xs, ys, pch = point_pch,
                     col = point_col_alpha,
                     cex = point_size, lwd = 1.2)
    if (isTRUE(label)) {
      graphics::text(xs, ys, labels = nm,
                     pos = 3, offset = 0.6,
                     cex = label_cex, col = "grey15")
    }

    if (p == 1L && (!is.null(grp_lev) || !is.null(shp_lev))) {
      leg_text <- c(grp_lev, shp_lev)
      leg_col  <- c(grp_cols[seq_along(grp_lev)],
                    rep("grey30", length(shp_lev)))
      leg_pch  <- c(rep(16L, length(grp_lev)),
                    shp_pchs[seq_along(shp_lev)])
      leg_title <- if (!is.null(grp_lev) && !is.null(shp_lev)) {
        "groups / shapes"
      } else if (!is.null(grp_lev)) {
        "groups"
      } else {
        "shapes"
      }
      graphics::legend(
        "topright", legend = leg_text,
        col = leg_col, pch = leg_pch,
        title = leg_title, title.col = "grey20",
        bty = "n", cex = 0.78, text.col = "grey20",
        inset = c(0.01, 0.01), pt.cex = 1.1
      )
    }
  }

  if (n_panels == 1L) {
    if (!is.null(main)) {
      graphics::title(main = main, line = 2.1,
                      cex.main = 1.0, col.main = "grey10")
    }
    if (isTRUE(show_gof)) {
      var_str <- paste(sprintf("%.1f%%", variance_pct),
                       collapse = " / ")
      gof <- sprintf("k = %d (%s)  stress = %.3f %s  var: %s",
                     k_use,
                     if (is_auto) "auto" else "user",
                     stress_1, band, var_str)
      graphics::mtext(gof, side = 3, line = 0.6,
                      cex = 0.78, col = "grey40")
    }
  } else {
    if (!is.null(main)) {
      graphics::mtext(main, side = 3, line = 1.6, outer = TRUE,
                      cex = 1.0, col = "grey10")
    }
    if (isTRUE(show_gof)) {
      var_str <- paste(sprintf("%.1f%%", variance_pct),
                       collapse = " / ")
      gof <- sprintf("k = %d (%s)  stress = %.3f %s  var: %s",
                     k_use,
                     if (is_auto) "auto" else "user",
                     stress_1, band, var_str)
      graphics::mtext(gof, side = 3, line = 0.2, outer = TRUE,
                      cex = 0.82, col = "grey40")
    }
  }

  out <- structure(
    list(
      mds_points         = mds_points,
      distance_matrix    = d_used,
      k_selected         = k_use,
      k_selection_mode   = if (is_auto) "auto" else "user",
      stress_by_k        = stress_by_k,
      variance_pct_by_k  = variance_pct_by_k,
      stress_threshold   = stress_threshold,
      stress_band        = band,
      eigenvalues        = mds$eig,
      variance_pct       = variance_pct,
      stress_1           = stress_1,
      distance           = distance,
      n_pixels_used      = n_pix_used,
      mask               = mask,
      groups             = groups,
      shapes             = shapes,
      panel_pairs        = panel_pairs,
      file               = file
    ),
    class = "rcisignal_mds"
  )
  invisible(out)
}

#' @keywords internal
#' @noRd
validate_categorical_arg <- function(x, ci_names, arg_name) {
  if (is.null(x)) return(NULL)
  if (!is.character(x) && !is.factor(x)) {
    cli::cli_abort(
      "{.arg {arg_name}} must be a named character vector."
    )
  }
  x_names <- names(x)
  x_vals  <- as.character(x)
  names(x_vals) <- x_names
  if (length(x_vals) != length(ci_names)) {
    cli::cli_abort(c(
      "{.arg {arg_name}} must have one entry per CI.",
      "i" = "Got length {length(x_vals)} vs {length(ci_names)} CIs."
    ))
  }
  if (is.null(x_names)) {
    cli::cli_abort(
      "{.arg {arg_name}} must be a NAMED character vector."
    )
  }
  if (!setequal(x_names, ci_names)) {
    missing <- setdiff(ci_names, x_names)
    extra   <- setdiff(x_names, ci_names)
    cli::cli_abort(c(
      "Names of {.arg {arg_name}} must match {.code names(cis)}.",
      if (length(missing)) c("x" = "Missing: {.val {missing}}.") else NULL,
      if (length(extra))   c("x" = "Extra: {.val {extra}}.") else NULL
    ))
  }
  unname(x_vals[ci_names])
}

#' @export
print.rcisignal_mds <- function(x, ...) {
  cat("<rcisignal MDS>\n")
  cat(sprintf("  n_cis:                %d\n", nrow(x$mds_points)))
  cat(sprintf("  mask:                 %s   (n_pixels_used: %s)\n",
              x$mask, format(x$n_pixels_used, big.mark = ",")))
  cat(sprintf("  distance:             %s\n", x$distance))
  cat("\n")
  cat(sprintf("  Dimensionality selection (%s, threshold = %.3f):\n",
              x$k_selection_mode, x$stress_threshold))
  k_grid <- as.integer(names(x$stress_by_k))
  for (i in seq_along(k_grid)) {
    s <- x$stress_by_k[[i]]
    v <- x$variance_pct_by_k[[i]]
    bnd <- if (s <= 0.025) "excellent"
           else if (s <= 0.05)  "good"
           else if (s <= 0.10)  "fair"
           else if (s <= 0.20)  "poor"
           else                 "very poor"
    arrow <- if (k_grid[i] == x$k_selected) " <--" else "    "
    cat(sprintf("    k = %d  stress = %.3f  %-9s%s   cum.var = %5.1f%%\n",
                k_grid[i], s, bnd, arrow, v))
  }
  cat("\n")
  cat(sprintf("  Selected: k = %d  (%s)\n",
              x$k_selected,
              if (x$k_selection_mode == "auto") {
                if (x$stress_1 <= x$stress_threshold) {
                  sprintf("smallest k reaching \"%s\"", x$stress_band)
                } else {
                  sprintf("fallback at k_max; best achievable stress = %.3f (%s)",
                          x$stress_1, x$stress_band)
                }
              } else {
                "user-specified"
              }))
  var_per_axis <- paste(
    sprintf("Dim %d %.1f%%", seq_along(x$variance_pct), x$variance_pct),
    collapse = "   "
  )
  cat(sprintf("  Variance per axis:    %s\n", var_per_axis))
  cat("\n")
  cat(sprintf("  Use `out$mds_points` to extract the %d x %d coordinate matrix\n",
              nrow(x$mds_points), ncol(x$mds_points)))
  cat("  in the Euclidean MDS space.\n")
  invisible(x)
}

#' @export
summary.rcisignal_mds <- function(object, ...) {
  print(object, ...)
  invisible(object)
}
