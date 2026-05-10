# Detect response-code inversion

Some participants (or whole batches of exported data) may have their
response codes flipped relative to the convention `rcicr` expects: `+1`
means "inverted version chosen" when the analyst intended "original
chosen", or the CSV export inadvertently negated the response column. A
straightforward way to detect this is to compute the infoVal twice per
participant — once with the original codes and once with every response
negated — and compare. For correctly coded data, the original should
score higher; if the flipped CI scores meaningfully higher, the response
column is probably inverted.

## Usage

``` r
check_response_inversion(
  responses,
  method = c("2ifc", "briefrc"),
  rdata = NULL,
  stimuli = NULL,
  baseimage = "base",
  col_participant = "participant_id",
  col_stimulus = "stimulus",
  col_response = "response",
  margin = 1.96,
  iter = 10000L,
  ...
)
```

## Arguments

- responses:

  Data frame with one row per trial. Required columns: `participant_id`,
  `stimulus`, `response` (values in `{-1, +1}`). Load yours from CSV via
  [`read_responses()`](https://olivethree.github.io/rcisignal/reference/read_responses.md)
  or [`utils::read.csv()`](https://rdrr.io/r/utils/read.table.html);
  column names are configurable via the `col_*` arguments.

- method:

  `"2ifc"` or `"briefrc"`.

- rdata:

  Path to the rcicr `.RData` file. Either `rdata` or `stimuli` must be
  supplied for the 2IFC path.

- stimuli:

  In-memory stimuli list (the `$stimuli` element of an `rcisignal_sim`
  object). Use in place of `rdata` when the file path no longer resolves
  (e.g. after [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)/
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html) across R sessions).

- baseimage:

  Name of the base image in `rdata$base_face_files`.

- col_participant, col_stimulus, col_response:

  Column names.

- margin:

  Numeric. Flagging threshold: flipped infoVal must exceed original
  infoVal by at least this amount. Default `1.96`.

- iter:

  Reference-distribution iterations. Default `10000`.

- ...:

  Passed through to
  [`compute_infoval_summary()`](https://olivethree.github.io/rcisignal/reference/compute_infoval_summary.md).

## Value

An
[`rcisignal_diag_result()`](https://olivethree.github.io/rcisignal/reference/rcisignal_diag_result.md)
object. `data$per_participant` has `participant_id`, `infoval_original`,
`infoval_flipped`, and `likely_inverted` (logical).

## Details

This check calls
[`compute_infoval_summary()`](https://olivethree.github.io/rcisignal/reference/compute_infoval_summary.md)
twice and reports participants whose flipped infoVal exceeds the
original by `margin` or more. A non-zero count is a strong signal that
the response column is miscoded for those participants.

Requires `rcicr` (same as
[`compute_infoval_summary()`](https://olivethree.github.io/rcisignal/reference/compute_infoval_summary.md)).
Runs two infoVal sweeps, so it takes roughly twice as long.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_2ifc_data(n_per_condition = 10, n_trials = 60, seed = 1)
check_response_inversion(sim$data, method = "2ifc",
                         rdata = sim$rdata_path)
} # }
```
