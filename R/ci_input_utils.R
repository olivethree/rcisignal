## Internal helpers shared by plot_ci_correlogram(),
## plot_ci_distance_matrix(), and plot_ci_mds(). Validate the
## `cis` argument (named list of CIs in any of the three accepted
## shapes), reduce per-producer matrices to group means, infer
## `img_dims`, and resolve a `mask` keyword to a logical vector.

#' @keywords internal
#' @noRd
prepare_ci_matrix <- function(cis, img_dims = NULL,
                              min_cis = 2L) {
  if (is.matrix(cis) && is.numeric(cis) &&
        !is.null(colnames(cis))) {
    src_dims <- attr(cis, "img_dims")
    cis <- stats::setNames(
      lapply(seq_len(ncol(cis)), function(j) cis[, j]),
      colnames(cis)
    )
    if (!is.null(src_dims) && is.null(img_dims)) {
      img_dims <- as.integer(src_dims)
    }
  }
  if (!is.list(cis)) {
    cli::cli_abort(c(
      "{.arg cis} must be a named list of CIs or a numeric matrix \\
       with named columns (one column per CI).",
      "i" = "Got: {.cls {class(cis)[1L]}}."
    ))
  }
  if (length(cis) < min_cis) {
    cli::cli_abort(
      "{.arg cis} must contain at least {min_cis} CIs."
    )
  }
  nm <- names(cis)
  if (is.null(nm) || any(!nzchar(nm))) {
    cli::cli_abort(c(
      "{.arg cis} must be fully named.",
      "i" = "Use {.code list(Trust = ..., Friendly = ...)} syntax."
    ))
  }

  reduce_one <- function(x, idx) {
    if (is.matrix(x) && is.numeric(x)) {
      if (ncol(x) == 1L) return(as.numeric(x))
      return(rowMeans(x))
    }
    if (is.numeric(x) && is.null(dim(x))) {
      return(as.numeric(x))
    }
    cli::cli_abort(
      "{.arg cis[[{idx}]]} must be a numeric vector or matrix."
    )
  }
  vecs <- lapply(seq_along(cis), function(i) reduce_one(cis[[i]], i))
  lens <- vapply(vecs, length, integer(1L))
  if (length(unique(lens)) != 1L) {
    cli::cli_abort(c(
      "All CIs in {.arg cis} must have the same number of pixels.",
      "i" = "Got lengths {.val {lens}}."
    ))
  }
  n_pix <- lens[1L]

  if (is.null(img_dims)) {
    a <- attr(cis[[1L]], "img_dims")
    if (!is.null(a)) {
      img_dims <- as.integer(a)
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
  if (length(img_dims) != 2L || any(img_dims < 1L)) {
    cli::cli_abort(
      "{.arg img_dims} must be a positive length-2 integer."
    )
  }
  if (prod(img_dims) != n_pix) {
    cli::cli_abort(c(
      "{.arg img_dims} disagrees with the CI length.",
      "i" = "prod(img_dims) = {prod(img_dims)}; CI length = {n_pix}."
    ))
  }

  M <- do.call(cbind, vecs)
  colnames(M) <- nm

  list(M = M, img_dims = img_dims, n_pix = n_pix, names = nm)
}

#' @keywords internal
#' @noRd
resolve_ci_mask <- function(mask, img_dims, n_pix) {
  if (identical(mask, "none")) {
    return(rep(TRUE, n_pix))
  }
  region <- if (identical(mask, "face")) "full" else mask
  as.logical(make_face_mask(img_dims, region = region))
}
