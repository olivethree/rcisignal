#' Read reverse-correlation response data from a CSV file
#'
#' @description
#' Reads a trial-level response CSV and returns a
#' [data.table::data.table]. Validates that the required columns are
#' present but does not coerce types or rename columns. Use this to
#' load responses for the diagnostics battery
#' ([run_diagnostics()] / `check_*()`) or as input to
#' [ci_from_responses_2ifc()] / [ci_from_responses_briefrc()].
#'
#' @param path Path to a CSV file. Any format readable by
#'   [data.table::fread()] is accepted (comma-, tab-, or
#'   space-delimited).
#' @param method Either `"2ifc"` or `"briefrc"`. Currently
#'   informational only; downstream checks use this to select
#'   method-specific logic.
#' @param col_participant,col_stimulus,col_response Column names that
#'   must exist in the file. Default to the package-wide standard
#'   names.
#' @param col_rt Optional column name for response time. If
#'   provided, must exist in the file.
#' @return A [data.table::data.table] with all columns from the CSV.
#' @seealso [run_diagnostics()], [ci_from_responses_2ifc()],
#'   [ci_from_responses_briefrc()]
#' @export
#' @examples
#' \dontrun{
#' responses <- read_responses("my_data.csv", method = "2ifc")
#' }
read_responses <- function(path,
                           method          = c("2ifc", "briefrc"),
                           col_participant = "participant_id",
                           col_stimulus    = "stimulus",
                           col_response    = "response",
                           col_rt          = NULL) {
  method <- match.arg(method)
  validate_path(path, "path")
  responses <- data.table::fread(path)
  validate_responses_df(
    responses,
    col_participant = col_participant,
    col_stimulus    = col_stimulus,
    col_response    = col_response,
    col_rt          = col_rt
  )
  responses
}
