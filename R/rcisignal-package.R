#' rcisignal: Quality Checks for Reverse-Correlation Data and
#' Classification Images
#'
#' @description
#' A consolidated toolkit for assessing the quality of reverse
#' correlation (RC) experiments end-to-end. Input-side diagnostics
#' (`run_diagnostics()` and the `check_*` family) catch silent
#' data-processing errors before classification image (CI)
#' computation. Output-side analyses (`run_reliability()`,
#' `run_discriminability()`, `infoval()`, `agreement_map_test()`)
#' quantify CI quality after computation.
#'
#' Works with both the standard 2IFC pipeline (via the upstream
#' `rcicr` package) and the Brief-RC pipeline (Schmitz, Rougier &
#' Yzerbyt, 2024 — implemented natively here). Operates directly on
#' the pixel-level signal produced by the original producers, so it
#' does not depend on a second-phase trait-rating study.
#'
#' @keywords internal
"_PACKAGE"
