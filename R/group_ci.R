#' Collapse a per-producer signal matrix into per-group means
#'
#' @description
#' Stage-2 aggregator of the rcisignal pipeline. Collapses a
#' per-producer `signal_matrix` (pixels x n_producers, the
#' object returned by [ci_from_responses_briefrc()] /
#' [ci_from_responses_2ifc()]) into a per-group matrix
#' (pixels x n_groups) for use with the distance-matrix, MDS, and
#' correlogram plot functions.
#'
#' @details
#' The package has two stages. Stage 1 (per-producer
#' `signal_matrix`) is the only object accepted by reliability,
#' discriminability, and infoVal functions. Stage 2
#' (group-averaged matrix) is for plotting, RDM-style comparison,
#' and MDS. `group_ci()` is the stage-1-to-stage-2 transformer.
#'
#' For each group, the corresponding output column is
#' `rowMeans(signal_matrix[, producers_in_group, drop = FALSE])`.
#' When `by` is a list, the cell label is the levels joined by
#' `"_"` in the order the list elements are given.
#'
#' `group_ci()` does not accept (and will never accept)
#' `trial_counts`, `noise_matrix`, or `mask`. Anything that needs
#' producer-level information lives in stage 1: do the analysis
#' first, then aggregate.
#'
#' @param signal_matrix Numeric matrix of pixels x n_producers, as
#'   returned by [ci_from_responses_briefrc()] or
#'   [ci_from_responses_2ifc()] in their `$signal_matrix` field.
#'   The column count must match `length(by)` (or
#'   `length(by[[1]])` when `by` is a list).
#' @param by Either an atomic vector of length
#'   `ncol(signal_matrix)` (one group label per producer) or a
#'   named list of such vectors for factorial grouping (cell
#'   names are the levels joined by `"_"`). Coerced to factor;
#'   `NA` levels are dropped with a warning naming the count.
#' @param drop Logical. When `TRUE` (default), empty cells are
#'   dropped from the output. When `FALSE`, empty cells are
#'   present as `NA` columns with `n = 0L`.
#'
#' @return A numeric matrix of pixels x n_groups, classed
#'   `c("rcisignal_group_ci", "matrix", "array")`. Carries:
#'   - column names = group labels;
#'   - `attr(., "n")` = named integer vector of per-group
#'     producer counts;
#'   - `attr(., "img_dims")` = inherited from
#'     `attr(signal_matrix, "img_dims")` if present.
#'
#' @seealso
#' Upstream (stage 1): [ci_from_responses_briefrc()],
#'   [ci_from_responses_2ifc()].
#' Downstream (stage 2): [plot_ci_distance_matrix()],
#'   [plot_ci_mds()], [plot_ci_correlogram()].
#'
#' @examples
#' set.seed(1)
#' n_pix  <- 32L * 32L
#' n_prod <- 12L
#' sm <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod,
#'              dimnames = list(NULL, sprintf("p%02d", seq_len(n_prod))))
#' attr(sm, "img_dims") <- c(32L, 32L)
#' g <- rep(c("A", "B"), each = n_prod / 2L)
#' gcis <- group_ci(sm, by = g)
#' dim(gcis)                 # n_pix x 2
#' attr(gcis, "n")           # named per-group producer counts
#'
#' \dontrun{
#' # Realistic end-to-end: simulate, build per-producer CIs, collapse.
#' sim <- simulate_briefrc_data(
#'   n_per_condition = 10, n_trials = 60,
#'   conditions = c("A", "B"), seed = 1
#' )
#' cis <- ci_from_responses_briefrc(sim$data,
#'                                  noise_matrix = sim$noise_matrix)
#' producer_to_cond <- sim$data$condition[match(
#'   colnames(cis$signal_matrix), sim$data$participant_id)]
#' gcis <- group_ci(cis$signal_matrix, by = producer_to_cond)
#' plot_ci_distance_matrix(gcis,
#'   img_dims = c(sim$meta$img_size, sim$meta$img_size))
#' }
#' @export
group_ci <- function(signal_matrix, by, drop = TRUE) {
  if (inherits(signal_matrix, "rcisignal_group_ci")) {
    cli::cli_abort(c(
      "{.arg signal_matrix} is already a stage-2 \\
       {.cls rcisignal_group_ci} object.",
      "i" = "{.fn group_ci} consumes a per-producer matrix \\
             (stage 1) and produces a group-averaged matrix \\
             (stage 2). Calling it on its own output would \\
             collapse twice."
    ))
  }
  if (!is.matrix(signal_matrix) || !is.numeric(signal_matrix)) {
    cli::cli_abort(
      "{.arg signal_matrix} must be a numeric matrix \\
       (pixels x n_producers)."
    )
  }
  n_prod <- ncol(signal_matrix)
  if (n_prod < 1L) {
    cli::cli_abort("{.arg signal_matrix} has no columns to group.")
  }

  factor_label <- build_grouping_factor(by, n_prod)
  f       <- factor_label$factor
  by_name <- factor_label$arg_text

  n_na <- sum(is.na(f))
  if (n_na > 0L) {
    cli::cli_warn(
      "Dropping {n_na} producer{?s} with NA in {.arg by}."
    )
  }

  groups <- levels(f)
  n_per  <- as.integer(tabulate(f, nbins = length(groups)))
  names(n_per) <- groups

  out <- matrix(NA_real_, nrow = nrow(signal_matrix),
                ncol = length(groups),
                dimnames = list(NULL, groups))
  for (k in seq_along(groups)) {
    idx <- which(as.integer(f) == k)
    if (length(idx) == 0L) next
    out[, k] <- rowMeans(signal_matrix[, idx, drop = FALSE])
  }

  if (isTRUE(drop)) {
    keep <- n_per > 0L
    if (any(!keep)) {
      out   <- out[, keep, drop = FALSE]
      n_per <- n_per[keep]
    }
  }

  attr(out, "n") <- n_per
  src_dims <- attr(signal_matrix, "img_dims")
  if (!is.null(src_dims)) attr(out, "img_dims") <- src_dims
  attr(out, "by_name") <- by_name
  class(out) <- c("rcisignal_group_ci", "matrix", "array")
  out
}

#' Build a single factor from an atomic vector or a list of vectors
#'
#' Returns a list with `factor` (the result) and `arg_text` (the
#' description used in print methods).
#'
#' @keywords internal
#' @noRd
build_grouping_factor <- function(by, n_prod) {
  if (is.list(by) && !is.data.frame(by)) {
    if (length(by) == 0L) {
      cli::cli_abort("{.arg by} list is empty.")
    }
    lens <- vapply(by, length, integer(1L))
    if (!all(lens == n_prod)) {
      cli::cli_abort(c(
        "{.arg by} list elements must each have length \\
         {n_prod} (= {.code ncol(signal_matrix)}).",
        "*" = "Got lengths: {.val {lens}}"
      ))
    }
    parts <- lapply(by, function(x) as.character(x))
    combined <- do.call(paste, c(parts, list(sep = "_")))
    has_na   <- Reduce(`|`, lapply(by, is.na))
    combined[has_na] <- NA_character_
    f <- factor(combined)
    if (!is.null(names(by))) {
      arg_text <- paste(names(by), collapse = " x ")
    } else {
      arg_text <- sprintf("list of %d factors", length(by))
    }
    return(list(factor = f, arg_text = arg_text))
  }
  if (length(by) != n_prod) {
    cli::cli_abort(c(
      "{.arg by} has length {length(by)} but \\
       {.code ncol(signal_matrix)} is {n_prod}.",
      "i" = "Pass one group label per producer."
    ))
  }
  f <- if (is.factor(by)) by else factor(by)
  list(factor = f, arg_text = "by")
}

#' @export
print.rcisignal_group_ci <- function(x, ...) {
  n_per <- attr(x, "n")
  n_grp <- ncol(x)
  n_pix <- nrow(x)
  by_text <- attr(x, "by_name") %||% "by"
  cat(sprintf("<rcisignal_group_ci: %d group%s, %d pixels>\n",
              n_grp, if (n_grp == 1L) "" else "s", n_pix))
  cat(sprintf("  grouping     : %s\n", by_text))
  if (!is.null(n_per)) {
    cat("  group sizes  :\n")
    show <- if (length(n_per) > 10L) c(head(n_per, 8L),
                                       stats::setNames(NA, "...")) else n_per
    for (nm in names(show)) {
      v <- show[[nm]]
      if (is.na(v)) {
        cat(sprintf("    %s (and %d more)\n",
                    nm, length(n_per) - 8L))
      } else {
        cat(sprintf("    %-20s n = %d\n", nm, v))
      }
    }
  }
  dims <- attr(x, "img_dims")
  if (!is.null(dims)) {
    cat(sprintf("  image dims   : %d x %d px\n", dims[1L], dims[2L]))
  }
  cat("Note: per-producer information has been averaged out.\n")
  cat("For reliability, infoVal, or discriminability, use the\n")
  cat("per-producer signal_matrix returned by ci_from_responses_*()\n")
  cat("(stage 1).\n")
  invisible(x)
}

#' @export
as.list.rcisignal_group_ci <- function(x, ...) {
  out <- lapply(seq_len(ncol(x)), function(j) x[, j])
  names(out) <- colnames(x)
  out
}

#' Abort if `x` is a stage-2 (group-averaged) CI matrix
#'
#' Two-line guard used at the top of every stage-1 function. The
#' error message names the calling function and points at the
#' two-stage pattern, so users who arrive via Stack Overflow /
#' Slack / vibes coding learn the architecture from the error
#' itself.
#'
#' @keywords internal
#' @noRd
abort_if_group_ci <- function(x, fn, arg = "signal_matrix") {
  if (!inherits(x, "rcisignal_group_ci")) {
    return(invisible(NULL))
  }
  cli::cli_abort(c(
    "{.fn {fn}} requires a per-producer CI matrix \\
     (stage 1), not a group-averaged matrix (stage 2).",
    "i" = "Pass the {.var signal_matrix} returned by \\
           {.fn ci_from_responses_briefrc} / \\
           {.fn ci_from_responses_2ifc} (one column per producer), \\
           not the matrix returned by {.fn group_ci} (one column \\
           per group).",
    "i" = "See {.code vignette(\"rcisignal\")} for the \\
           two-stage pattern."
  ),
  call = NULL)
}
