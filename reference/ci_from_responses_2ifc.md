# Compute individual 2IFC CIs from trial-level responses

Thin wrapper around
[`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html)
that handles the known rcicr-integration gotchas and exposes a uniform
return shape across the 2IFC and Brief-RC paths.

Use this when you have 2IFC trial-level responses and want the package
to compute individual CIs ready for the reliability metrics.

## Usage

``` r
ci_from_responses_2ifc(
  responses,
  rdata_path = NULL,
  stimuli = NULL,
  base_image = NULL,
  col_participant = "participant_id",
  col_stimulus = "stimulus",
  col_response = "response",
  scaling = "autoscale",
  keep_rendered = FALSE,
  targetpath = tempfile("rcisignal_2ifc_"),
  save_as_png = FALSE
)
```

## Arguments

- responses:

  A data frame with one row per trial. Must contain `col_participant`
  (producer id), the column named by `col_stimulus` (stimulus id,
  integer, index into the rcicr noise pool), and `col_response` with
  values in `{-1, +1}`.

- rdata_path:

  Path to the `.Rdata` file produced by
  [`rcicr::generateStimuli2IFC()`](https://rdrr.io/pkg/rcicr/man/generateStimuli2IFC.html).
  Either `rdata_path` or `stimuli` must be supplied; if both are given
  `stimuli` wins.

- stimuli:

  In-memory stimuli list as returned in `$stimuli` by
  [`simulate_2ifc_data()`](https://olivethree.github.io/rcisignal/reference/simulate_2ifc_data.md).
  Use this in place of `rdata_path` when the simulation has been saved
  with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and reloaded
  in a different R session: the path stored on `$rdata_path` no longer
  resolves, but `$stimuli` is self-contained. Internally the list is
  written to a fresh tempdir-backed `.Rdata` before the call into rcicr.

- base_image:

  Base face image. Three forms are accepted, matching
  [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md):

  - a numeric matrix in `[0, 1]` (used directly),

  - a single string path to a PNG / JPEG, or

  - a single string label naming an entry in the rdata's `base_faces`
    list (the historical 2IFC form). When `NULL`, the rdata's single
    base is used; if the rdata contains more than one base, the call
    aborts and lists the available labels. Matrix and path forms are
    injected into a temporary copy of the rdata under a synthetic label
    so the rcicr call sees the same structure it always has.

- col_participant, col_stimulus, col_response:

  Column names in `responses`.

- scaling:

  rcicr scaling option; one of `"autoscale"`, `"independent"`,
  `"constant"`, `"none"`. Passed through to
  [`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html).
  Scaling only affects the `$combined` (rendered) image, not `$ci`, so
  the returned `$signal_matrix` is the raw mask regardless of this
  argument.

- keep_rendered:

  If `TRUE`, also extract rcicr's `$combined` image (base + scaled
  noise) per producer and stack into `$rendered_ci`. Default `FALSE`.
  Visualization only.

- targetpath:

  Where rcicr writes PNGs. Defaults to an auto-cleaned tempdir so the
  working directory is not polluted.

- save_as_png:

  Whether rcicr writes individual CI PNGs. Defaults to `FALSE` for speed
  in a pure reliability pipeline.

## Value

A list with `signal_matrix`, optionally `rendered_ci`, `participants`,
`img_dims`, `scaling`, and `rcicr_result`.

## Details

What the wrapper does for you:

- Attaches `foreach`, `tibble`, `dplyr` at runtime (`%dopar%` /
  `tribble` / `%>%` are not namespace-prefixed inside rcicr).

- Matches the `.Rdata` extension case-insensitively (rcicr writes
  lowercase on some filesystems).

- Extracts the per-participant CI noise component (`$ci`) from the rcicr
  result list and stacks it into a pixels x participants signal matrix,
  already base-subtracted, ready for `rel_*()`.

- When `keep_rendered = TRUE`, also extracts the rendered `$combined`
  image and stacks it as `$rendered_ci` for visualization (display only,
  not for downstream stats).

`rcicr` must be installed (it is a Suggests dep; install with
`remotes::install_github("rdotsch/rcicr")` if needed).

## Reading the result

- `$signal_matrix` is the **raw mask** (rcicr's `$ci` per producer).
  This is what every `rel_*` function expects.

- `$rendered_ci`, present when `keep_rendered = TRUE`, is the
  `base + scaling(mask)` image rcicr would have written to PNG.
  **Visualization only.**

- `$rcicr_result` is the raw return value of
  [`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html).

## Common mistakes

- Passing `$rendered_ci` to `rel_*` functions; it carries the scaling
  step, which distorts variance-based metrics.

- Forgetting that the warning about Mode 1 also applies if you later
  read `$rendered_ci` back from PNG instead of using `$signal_matrix`
  directly.

## See also

[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md),
[`rcicr::batchGenerateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/batchGenerateCI2IFC.html),
[`group_ci()`](https://olivethree.github.io/rcisignal/reference/group_ci.md)
(stage 2, optional)

## Examples

``` r
if (FALSE) { # \dontrun{
# In your study, you replace the simulated inputs with your own:
#   - `sim$data` -> your responses data frame (CSV via read_responses())
#   - `sim$rdata_path` -> path to the rcicr stimuli .Rdata file used to
#     generate your stimuli.
# The `simulate_2ifc_data()` helper just produces the same shapes so the
# example below is runnable end-to-end.
sim <- simulate_2ifc_data(n_per_condition = 10, n_trials = 60, seed = 1)
res <- ci_from_responses_2ifc(sim$data, rdata_path = sim$rdata_path)
dim(res$signal_matrix)   # n_pixels x n_producers
rel_split_half(res$signal_matrix, n_permutations = 200L, seed = 1)

# Same call, but using the portable in-memory stimuli list -- the
# form that survives saveRDS()/readRDS() across R sessions.
res2 <- ci_from_responses_2ifc(sim$data, stimuli = sim$stimuli)

# Supply a base image as a path or matrix (same shapes as Brief-RC).
res3 <- ci_from_responses_2ifc(sim$data, rdata_path = sim$rdata_path,
                               base_image = sim$base_image_path)
res4 <- ci_from_responses_2ifc(sim$data, rdata_path = sim$rdata_path,
                               base_image = sim$base_face)
} # }
```
