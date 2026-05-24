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
#' Uses the same "data frame plus column name" idiom as every
#' other responses-consuming function in the package: pass the
#' trial-level `responses` data frame, plus the name of the
#' column you want to group by. Producer-to-group alignment
#' happens internally via `colnames(signal_matrix)`.
#'
#' @details
#' The package has two stages. Stage 1 (per-producer
#' `signal_matrix`) is the only object accepted by reliability,
#' discriminability, and informational-value functions. Stage 2
#' (group-averaged matrix) is for plotting, RDM-style comparison,
#' and MDS. `group_ci()` is the stage-1-to-stage-2 transformer.
#'
#' For each group, the corresponding output column is
#' `rowMeans(signal_matrix[, producers_in_group, drop = FALSE])`.
#' For a factorial `by` (length 2+), the cell label is the levels
#' joined by `"_"` in the column order given.
#'
#' `group_ci()` does not accept (and will never accept)
#' `trial_counts`, `noise_matrix`, or `mask`. Anything that needs
#' producer-level information lives in stage 1: do the analysis
#' first, then aggregate.
#'
#' @param signal_matrix Numeric matrix of pixels x n_producers, as
#'   returned by [ci_from_responses_briefrc()] or
#'   [ci_from_responses_2ifc()] in their `$signal_matrix` field.
#'   Must have non-empty column names (the producer ids that map
#'   into `responses[[col_participant]]`).
#' @param responses Trial-level data frame containing one row per
#'   trial. Must contain the column named by `col_participant`
#'   (mapping to `colnames(signal_matrix)`) and every column
#'   named in `by`. Each producer's `by` value(s) must be
#'   consistent across all of their rows; an inconsistency aborts
#'   with a teaching message naming the offending producer.
#' @param by Character vector of column names in `responses`.
#'   Length 1 selects a single grouping column (e.g.
#'   `by = "condition"`). Length 2+ produces a factorial grouping
#'   where cell labels are the levels joined with `"_"` in the
#'   given order (e.g. `by = c("condition", "sex")` yields cells
#'   like `"happy_F"`).
#' @param col_participant Name of the participant-id column in
#'   `responses`. Default `"participant_id"`, matching the rest
#'   of the package.
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
#' @aliases rcisignal_group_ci
#'
#' @examples
#' \dontrun{
#' # End-to-end: simulate two conditions, build per-producer CIs,
#' # collapse into one CI per condition.
#' sim <- simulate_briefrc_data(
#'   n_per_condition = 10, n_trials = 60,
#'   conditions = c("A", "B"), seed = 1
#' )
#' cis  <- ci_from_responses_briefrc(sim$data,
#'                                   noise_matrix = sim$noise_matrix)
#' gcis <- group_ci(cis$signal_matrix, sim$data, by = "condition")
#' gcis                              # n_pixels x 2 (A and B)
#' attr(gcis, "n")                   # per-group producer counts
#'
#' # Factorial grouping via two columns:
#' # gcis_fact <- group_ci(cis$signal_matrix, sim$data,
#' #                       by = c("condition", "sex"))
#' }
#' @export
group_ci <- function(signal_matrix,
                     responses,
                     by,
                     col_participant = "participant_id",
                     drop            = TRUE) {
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
  if (ncol(signal_matrix) < 1L) {
    cli::cli_abort("{.arg signal_matrix} has no columns to group.")
  }
  if (is.null(colnames(signal_matrix))) {
    cli::cli_abort(c(
      "{.arg signal_matrix} has no column names.",
      "i" = "Column names are the producer ids that map into \\
             {.code responses[[col_participant]]}. The \\
             {.fn ci_from_responses_briefrc} / \\
             {.fn ci_from_responses_2ifc} builders set these for \\
             you; if you built {.arg signal_matrix} by hand, set \\
             {.code colnames(signal_matrix)} to your producer ids \\
             before calling."
    ))
  }

  factor_label <- build_grouping_from_responses(
    signal_matrix, responses, by, col_participant
  )
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

#' Build the per-producer grouping factor from a responses data frame
#'
#' Returns a list with `factor` (length `ncol(signal_matrix)`,
#' one entry per producer, in colnames order) and `arg_text` (the
#' description used in print methods).
#'
#' @keywords internal
#' @noRd
build_grouping_from_responses <- function(signal_matrix, responses,
                                          by, col_participant) {
  if (!is.data.frame(responses)) {
    cli::cli_abort(
      "{.arg responses} must be a data frame (got \\
       {.cls {class(responses)[1L]}})."
    )
  }
  if (!is.character(by) || length(by) < 1L || any(is.na(by)) ||
        any(!nzchar(by))) {
    cli::cli_abort(c(
      "{.arg by} must be a character vector of column names \\
       in {.arg responses}.",
      "i" = "Single grouping: {.code by = \"condition\"}.",
      "i" = "Factorial: {.code by = c(\"condition\", \"sex\")}."
    ))
  }
  required_cols <- c(col_participant, by)
  missing_cols  <- setdiff(required_cols, names(responses))
  if (length(missing_cols) > 0L) {
    n_missing <- length(missing_cols)
    cli::cli_abort(c(
      "{n_missing} column{?s} not in {.arg responses}: \\
       {.val {missing_cols}}.",
      "i" = "Available: {.val {names(responses)}}"
    ))
  }

  producer_ids <- colnames(signal_matrix)
  resp_pids    <- as.character(responses[[col_participant]])
  missing_pids <- setdiff(producer_ids, resp_pids)
  if (length(missing_pids) > 0L) {
    n_show <- min(5L, length(missing_pids))
    cli::cli_abort(c(
      "{length(missing_pids)} producer{?s} in \\
       {.code colnames(signal_matrix)} not found in \\
       {.code responses[[\"{col_participant}\"]]}.",
      "*" = "Missing: {.val {missing_pids[seq_len(n_show)]}}\\
             {if (length(missing_pids) > n_show) ', ...' else ''}",
      "i" = "{.arg signal_matrix} columns and {.arg responses} \\
             rows must come from the same dataset."
    ))
  }

  group_strings <- character(length(producer_ids))
  for (i in seq_along(producer_ids)) {
    pid     <- producer_ids[i]
    pid_rows <- which(resp_pids == pid)
    if (length(by) == 1L) {
      vals <- as.character(responses[[by]][pid_rows])
    } else {
      vals_per_row <- lapply(by, function(b) {
        as.character(responses[[b]][pid_rows])
      })
      vals <- do.call(paste, c(vals_per_row, list(sep = "_")))
    }
    uniq <- unique(vals)
    if (length(uniq) > 1L) {
      cli::cli_abort(c(
        "Producer {.val {pid}} has inconsistent {.arg by} \\
         values across their rows in {.arg responses}.",
        "*" = "Got: {.val {uniq}}",
        "i" = "Each producer's {.arg by} column(s) must be \\
               constant across their trials. Check the column{?s} \\
               {.val {by}} in {.arg responses}."
      ))
    }
    group_strings[i] <- if (length(uniq) == 0L) NA_character_ else uniq
  }

  # Preserve empty factor levels if `by` is a single factor column;
  # this is what makes `drop = FALSE` produce NA columns for unused
  # levels. For factorial `by`, levels are observed only.
  if (length(by) == 1L && is.factor(responses[[by]])) {
    src_levels <- levels(responses[[by]])
    all_levels <- union(stats::na.omit(unique(group_strings)),
                        src_levels)
    f <- factor(group_strings, levels = all_levels)
  } else {
    f <- factor(group_strings)
  }
  arg_text <- if (length(by) == 1L) by else paste(by, collapse = " x ")
  list(factor = f, arg_text = arg_text)
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
