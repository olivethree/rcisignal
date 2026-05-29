# Detect response-code inversion

Some participants (or whole batches of exported data) may have their
response codes flipped relative to the convention `rcicr` expects: `+1`
means "inverted version chosen" when the analyst intended "original
chosen", or the CSV export inadvertently negated the response column. A
straightforward way to detect this is to compute the infoVal twice per
participant – once with the original codes and once with every response
negated – and compare. For correctly coded data, the original should
score higher; if the flipped CI scores meaningfully higher, the response
column is probably inverted.

## Usage

``` r
check_response_inversion(
  responses,
  method = c("2ifc", "briefrc"),
  rdata = NULL,
  stimuli = NULL,
  noise_matrix = NULL,
  base_image = "base",
  col_participant = "participant_id",
  col_stimulus = "stimulus",
  col_response = "response",
  margin = 1.96,
  iter = 1000L,
  seed = NULL,
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

  Path to the rcicr `.RData` file (2IFC). Either `rdata` or `stimuli`
  must be supplied for the 2IFC path.

- stimuli:

  In-memory stimuli list (the `$stimuli` element of an `rcisignal_sim`
  object). Use in place of `rdata` when the file path no longer resolves
  (e.g. after [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)/
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html) across R sessions).

- noise_matrix:

  Noise matrix for the Brief-RC path. Either a numeric matrix of
  `n_pixels x pool_size`, or a path to a text / `.rds` file (see
  [`read_noise_matrix()`](https://olivethree.github.io/rcisignal/reference/read_noise_matrix.md)).

- base_image:

  Name of the base image in `rdata$base_face_files` (2IFC only). Default
  `"base"`.

- col_participant, col_stimulus, col_response:

  Column names.

- margin:

  Numeric. Flagging threshold: flipped infoVal must exceed original
  infoVal by at least this amount. Default `1.96`.

- iter:

  Reference-distribution iterations. Default `1000L`.

- seed:

  Optional integer; RNG state restored on exit.

- ...:

  Unused.

## Value

An
[`rcisignal_diag_result()`](https://olivethree.github.io/rcisignal/reference/rcisignal_diag_result.md)
object. `data$per_participant` has `participant_id`, `infoval_original`,
`infoval_flipped`, and `likely_inverted` (logical).

## Details

This check computes per-producer infoVal z-scores via the package-native
[`infoval()`](https://olivethree.github.io/rcisignal/reference/infoval.md)
pipeline twice (original and sign-flipped) and reports participants
whose flipped infoVal exceeds the original by `margin` or more. Supports
both 2IFC and Brief-RC. A non-zero count is a strong signal that the
response column is miscoded for those participants. Runs two infoVal
sweeps, so it takes roughly twice as long.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_2ifc_data(n_per_condition = 10, n_trials = 60, seed = 1)
check_response_inversion(sim$data, method = "2ifc",
                         rdata = sim$rdata_path)
} # }
```
