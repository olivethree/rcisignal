# Read reverse-correlation response data from a CSV file

Reads a trial-level response CSV and returns a
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html).
Validates that the required columns are present but does not coerce
types or rename columns. Use this to load responses for the diagnostics
battery
([`run_diagnostics()`](https://olivethree.github.io/rcisignal/reference/run_diagnostics.md)
/ `check_*()`) or as input to
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md)
/
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md).

## Usage

``` r
read_responses(
  path,
  method = c("2ifc", "briefrc"),
  col_participant = "participant_id",
  col_stimulus = "stimulus",
  col_response = "response",
  col_rt = NULL
)
```

## Arguments

- path:

  Path to a CSV file. Any format readable by
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)
  is accepted (comma-, tab-, or space-delimited).

- method:

  Either `"2ifc"` or `"briefrc"`. Currently informational only;
  downstream checks use this to select method-specific logic.

- col_participant, col_stimulus, col_response:

  Column names that must exist in the file. Default to the package-wide
  standard names.

- col_rt:

  Optional column name for response time. If provided, must exist in the
  file.

## Value

A
[data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
with all columns from the CSV.

## See also

[`run_diagnostics()`](https://olivethree.github.io/rcisignal/reference/run_diagnostics.md),
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md),
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
responses <- read_responses("my_data.csv", method = "2ifc")
} # }
```
