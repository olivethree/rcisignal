# Cross-validate infoVal against RT quality

Correlates per-participant infoVal with per-participant median response
time. The check looks for two patterns that are individually plausible
but jointly suspicious:

## Usage

``` r
check_rt_infoval_consistency(
  responses,
  method = c("2ifc", "briefrc"),
  rdata = NULL,
  stimuli = NULL,
  baseimage = "base",
  col_participant = "participant_id",
  col_stimulus = "stimulus",
  col_response = "response",
  col_rt = "rt",
  iter = 10000L,
  ...
)
```

## Arguments

- responses:

  Data frame with one row per trial. Required columns: `participant_id`,
  `stimulus`, `response` (values in `{-1, +1}`), and `rt` (response
  time, milliseconds). Load yours from CSV via
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

- col_participant, col_stimulus, col_response, col_rt:

  Column names. `col_rt` is required.

- iter:

  Reference-distribution iterations. Default `10000`.

- ...:

  Unused.

## Value

An
[`rcisignal_diag_result()`](https://olivethree.github.io/rcisignal/reference/rcisignal_diag_result.md)
object. `data$per_participant` has `participant_id`, `infoval`,
`median_rt`. `data$correlation` is the Pearson correlation between the
two.

## Details

1.  High infoVal with fast median RT. A participant producing a
    seemingly informative CI while responding faster than others is more
    likely to have produced a spurious signal than genuine meaningful
    information — the signal may be an artefact of a button-mashing
    strategy that happens to correlate with the noise.

2.  A negative correlation between infoVal and median RT across
    participants. If fast responders systematically score higher on
    infoVal, something about the measurement is off.

The function computes `cor(infoval, median_rt)` across participants and
returns the value plus a per-participant table for further inspection.
It does not auto-exclude anyone; interpretation requires judgement about
the specific experiment.

Requires `rcicr` via
[`compute_infoval_summary()`](https://olivethree.github.io/rcisignal/reference/compute_infoval_summary.md).

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_2ifc_data(n_per_condition = 10, n_trials = 60, seed = 1)
check_rt_infoval_consistency(sim$data, method = "2ifc",
                             rdata = sim$rdata_path, col_rt = "rt")
} # }
```
