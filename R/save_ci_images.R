#' Write CIs to PNG or JPEG files
#'
#' @description
#' Writes each column of a `signal_matrix` (one CI) to disk as its own
#' image. The default output matches what `rcicr::generateCI()` /
#' `rcicr::generateCI2IFC()` would write for the same CI: a grayscale
#' luminance image of the CI noise scaled into `[0, 1]` and averaged
#' with the base face (no color palette). Two palette overlays are
#' available as opt-ins for visualization: `"diverging"` (matches
#' `plot_ci_overlay()`, signed signal, blue = positive, red = negative)
#' and `"fire"` (matches `plot_agreement_map(palette = "fire")`,
#' unipolar `|t|`-style yellow-to-red).
#'
#' Works for both per-producer matrices (one column per producer, the
#' `$signal_matrix` returned by `ci_from_responses_*()`) and
#' group-level matrices (one column per group, the output of
#' [group_ci()] or `$group_ci` when the generator is called with
#' `group_by =`). Filenames are derived from the column names of
#' `signal_matrix`.
#'
#' @details
#' Filenames default to `<prefix><colname>.<ext>`, where `<prefix>` is
#' chosen automatically from `attr(signal_matrix, "ci_level")`:
#'
#' * `"individual"` (set by `ci_from_responses_*()` on the per-producer
#'   `$signal_matrix`) -> `prefix = "ind_ci_"`.
#' * `"group"` (set by [group_ci()] on its return matrix) ->
#'   `prefix = "group_ci_"`.
#' * No `ci_level` attribute -> defaults to `prefix = "ind_ci_"`.
#'
#' Override the auto-prefix by passing `prefix =` explicitly.
#'
#' The default `palette = "grayscale"` reproduces rcicr's
#' `generateCI(..., save_as_png = TRUE)` output exactly: for each CI
#' column, the raw noise is scaled into `[0, 1]` via the chosen
#' `scaling` method (default `"independent"`, matching rcicr's
#' default), then averaged with the base via `(scaled + base) / 2`,
#' then written via `png::writePNG()` (or `jpeg::writeJPEG()`) as a
#' grayscale image. Pass `palette = "diverging"` or `"fire"` instead
#' to write a colored overlay rendered the same way as the on-screen
#' plot functions.
#'
#' @param signal_matrix Numeric matrix with non-empty, unique column
#'   names. Per-producer or group-level; both are accepted.
#' @param base_image Base face image. Either a numeric matrix in
#'   `[0, 1]` or a single string path to a PNG / JPEG. Used as the
#'   underlay for every rendered CI.
#' @param dir Output directory. Created (recursively) if missing.
#' @param format Output format. `"png"` (default) or `"jpeg"`.
#' @param palette Color palette. `"grayscale"` (default; raw pixel
#'   luminance, matches rcicr), `"diverging"` (signed CI on a
#'   blue/red ramp, matches `plot_ci_overlay()`), or `"fire"`
#'   (unipolar `|t|`-style yellow-to-red).
#' @param scaling Scaling method for the `"grayscale"` palette,
#'   matching rcicr's `generateCI(scaling = ...)`: `"independent"`
#'   (default; per-CI symmetric scaling by `max(|ci|)`), `"constant"`
#'   (scale by a user-supplied `scaling_constant`, comparable across
#'   CIs), `"matched"` (range-match each CI to the base image range),
#'   or `"none"` (write the raw `ci + base` with no scaling; rarely
#'   what you want). Ignored when `palette != "grayscale"`.
#' @param scaling_constant Numeric constant used when
#'   `scaling = "constant"`. Default `0.1`, matching rcicr.
#' @param prefix Optional character scalar overriding the
#'   auto-derived filename prefix (`"ind_ci_"` for per-producer
#'   matrices, `"group_ci_"` for group-level). Pass any string to
#'   force a custom convention (e.g. `prefix = "trust_"`).
#' @param threshold Optional numeric. Pixels with absolute CI value
#'   below `threshold` are forced to 0 (grayscale) or to neutral
#'   (palette overlays).
#' @param mask Optional logical vector of length `nrow(signal_matrix)`.
#'   Pixels with `mask = FALSE` are set to `NA` (grayscale, matching
#'   rcicr's `applyMask()` semantics) or rendered as base only
#'   (palette overlays).
#' @param zlim Optional `c(lo, hi)` color-scale endpoints. Used only
#'   for `palette = "diverging"` and `palette = "fire"`. Ignored for
#'   grayscale.
#' @param alpha_max Numeric in `[0, 1]`. Maximum opacity of the
#'   heatmap at the color-scale top. Used only for palette overlays.
#'   Default `0.7`.
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
#' # Default: rcicr-style grayscale (raw luminance, no palette).
#' # Files: ind_ci_P001.png, ind_ci_P002.png, ...
#' save_ci_images(res$signal_matrix, base_image = sim$base_face,
#'                dir = out)
#'
#' # Group-level CIs, same rcicr-style grayscale output.
#' save_ci_images(res$group_ci, base_image = sim$base_face,
#'                dir = out)
#'
#' # Diverging blue/red overlay (rcisignal visualization, not rcicr).
#' save_ci_images(res$group_ci, base_image = sim$base_face,
#'                dir = out, palette = "diverging",
#'                prefix = "diverging_")
#' }
#' @export
save_ci_images <- function(signal_matrix,
                           base_image,
                           dir,
                           format           = c("png", "jpeg"),
                           palette          = c("grayscale",
                                                "diverging", "fire"),
                           scaling          = c("independent",
                                                "constant",
                                                "matched", "none"),
                           scaling_constant = 0.1,
                           prefix           = NULL,
                           threshold        = NULL,
                           mask             = NULL,
                           zlim             = NULL,
                           alpha_max        = 0.7,
                           img_dims         = NULL,
                           quality          = 90,
                           overwrite        = FALSE,
                           quiet            = FALSE) {
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
  scaling <- match.arg(scaling)
  ext <- switch(format, png = "png", jpeg = "jpeg")

  pkg <- switch(format, png = "png", jpeg = "jpeg")
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg {pkg}} is required to write {format} files.",
      "i" = "Install with {.code install.packages(\"{pkg}\")}."
    ))
  }

  if (is.null(prefix)) {
    level <- attr(signal_matrix, "ci_level")
    prefix <- if (identical(level, "group")) "group_ci_" else "ind_ci_"
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

  base_vec <- as.vector(base_mat)

  for (j in seq_len(ncol(signal_matrix))) {
    sig <- signal_matrix[, j]

    if (palette == "grayscale") {
      ci_vec <- sig
      if (!is.null(threshold)) {
        ci_vec[abs(ci_vec) < threshold] <- 0
      }
      if (!is.null(mask)) {
        ci_vec[!mask] <- NA_real_
      }
      scaled <- apply_rcicr_scaling(base_vec, ci_vec,
                                    scaling, scaling_constant)
      combined <- (scaled + base_vec) / 2
      combined[combined < 0 & !is.na(combined)] <- 0
      combined[combined > 1 & !is.na(combined)] <- 1
      combined_mat <- matrix(combined, img_dims[1L], img_dims[2L])

      if (format == "png") {
        png::writePNG(combined_mat, target = paths[j])
      } else {
        jpeg::writeJPEG(combined_mat, target = paths[j],
                        quality = quality / 100)
      }
      next
    }

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

#' Apply rcicr-style scaling to a raw CI noise vector
#'
#' Ports `rcicr:::applyScaling()` so `save_ci_images(palette =
#' "grayscale")` produces the same scaled noise rcicr would write.
#' `base_vec` and `ci` are flat numeric vectors over the same pixel
#' grid; `ci` may contain `NA`s (from `mask = FALSE`).
#'
#' @keywords internal
#' @noRd
apply_rcicr_scaling <- function(base_vec, ci, scaling, constant) {
  switch(
    scaling,
    none = ci,
    constant = {
      out <- (ci + constant) / (2 * constant)
      finite_out <- out[is.finite(out)]
      if (length(finite_out) &&
            (max(finite_out) > 1 || min(finite_out) < 0)) {
        cli::cli_warn(c(
          "{.code scaling_constant = {constant}} produced pixels \\
           outside {.code [0, 1]}; clipping will occur.",
          "i" = "Pick a larger {.arg scaling_constant} to avoid \\
                 clipping (matches rcicr's behavior)."
        ))
      }
      out
    },
    matched = {
      base_min <- min(base_vec, na.rm = TRUE)
      base_max <- max(base_vec, na.rm = TRUE)
      ci_finite <- ci[!is.na(ci)]
      ci_min <- min(ci_finite)
      ci_max <- max(ci_finite)
      ci_rng <- ci_max - ci_min
      if (!is.finite(ci_rng) || ci_rng == 0) {
        rep(base_min, length(ci))
      } else {
        base_min + (base_max - base_min) * (ci - ci_min) / ci_rng
      }
    },
    independent = {
      ci_finite <- ci[!is.na(ci)]
      r <- range(ci_finite)
      k <- max(abs(r))
      if (!is.finite(k) || k == 0) k <- 1
      (ci + k) / (2 * k)
    }
  )
}
