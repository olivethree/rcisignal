## Per-result-class regression list. When `print.*()` is called on
## an rcisignal result whose `attr(., "rcisignal_version")` matches a
## known-bad version for that class, emit a one-time per-session
## warning telling the user to recompute.
##
## Add an entry by listing it directly in `.rcisignal_regressions`
## below; the existing entry is the rcicrely-pre-merge bug
## preserved so users with v0.2.0 result objects on disk still get
## warned when they print them under rcisignal.
##
## Each entry: list(class = <result class>,
##                  version_pkg = "rcicrely" | "rcisignal",
##                  version = <package_version>,
##                  message = <user-facing string>).

.rcisignal_regressions <- list(
  list(
    class       = "rcisignal_rel_infoval",
    version_pkg = "rcicrely",
    version     = package_version("0.2.0"),
    message     = "Result was computed by rcicrely 0.2.0, which \\
                   sampled reference stim ids with replacement. \\
                   This biased z-scores downward (typically by ~1 \\
                   unit at n_trials == n_pool). Recompute under \\
                   rcisignal for correct calibration."
  )
)

#' Lookup table of known regressions for a given result class
#'
#' Returns the list of `(version, message)` pairs registered for the
#' given S3 class, or `NULL` if none.
#'
#' @keywords internal
#' @noRd
known_regressions_for <- function(cls) {
  hits <- vapply(.rcisignal_regressions,
                 function(e) identical(e$class, cls),
                 logical(1L))
  if (!any(hits)) return(NULL)
  .rcisignal_regressions[hits]
}

#' Emit the once-per-session "result computed by buggy version" warning
#'
#' Called from `print.*()` methods after looking up the result's
#' `rcisignal_version` attribute against `.rcisignal_regressions`.
#'
#' @keywords internal
#' @noRd
warn_known_regression <- function(x) {
  if (isTRUE(getOption("rcisignal.silence_regression_warning",
                       FALSE))) {
    return(invisible())
  }
  cls <- class(x)[[1L]]
  v   <- attr(x, "rcisignal_version")
  if (is.null(v)) return(invisible())
  hits <- known_regressions_for(cls)
  if (is.null(hits)) return(invisible())
  for (h in hits) {
    if (identical(v, h$version)) {
      key <- paste0("regression_", cls, "_", as.character(v))
      emitted <- isTRUE(.rcisignal_session_state[[key]])
      if (emitted) next
      pkg <- if (is.null(h$version_pkg)) "rcisignal" else h$version_pkg
      cli::cli_warn(c(
        "{.cls {cls}} result was computed by {pkg} \\
         {as.character(v)} (known regression).",
        "i" = h$message,
        "i" = "Silence: \\
               {.code options(rcisignal.silence_regression_warning = TRUE)}."
      ))
      .rcisignal_session_state[[key]] <- TRUE
    }
  }
  invisible()
}
