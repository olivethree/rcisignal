# Compute individual Brief-RC CIs from trial-level responses

Native implementation of Schmitz, Rougier & Yzerbyt (2024)'s Brief-RC
mask. Does not call any `rcicr::*_brief` function; those do not exist in
upstream rcicr v1.0.1. Only rcicr's noise-pattern pool (the `stimuli`
object inside an `.Rdata` from `generateStimuli2IFC()`) is reused; the
mask is computed in pure R.

Use this when you have Brief-RC 12 or Brief-RC 20 trial-level responses
and want the package to produce per-producer noise masks ready for the
reliability metrics.

## Usage

``` r
ci_from_responses_briefrc(
  responses,
  rdata_path = NULL,
  noise_matrix = NULL,
  base_image = NULL,
  col_participant = "participant_id",
  col_stimulus = "stimulus",
  col_response = "response",
  method = c("briefrc12", "briefrc20"),
  scaling = c("none", "matched", "constant"),
  scaling_constant = NULL,
  base_image_path = NULL
)
```

## Arguments

- responses:

  Data frame with one row per trial. Must contain the columns named by
  `col_participant`, `col_stimulus`, `col_response`. `response` values
  must be in `{-1, +1}`.

- rdata_path, noise_matrix:

  Exactly one must be supplied. Provide `rdata_path` to read the noise
  matrix from an rcicr `.Rdata`, or pass a pre-loaded `noise_matrix`
  directly.

- base_image:

  Base face image. Either a numeric matrix in `[0, 1]` (e.g.
  `sim$base_face` from
  [`simulate_briefrc_data()`](https://olivethree.github.io/rcisignal/reference/simulate_briefrc_data.md),
  or a hand-built mask) or a single string path to a PNG / JPEG.
  Optional when `scaling = "none"` (default), in which case the image is
  not needed for the mathematical mask. Required for the
  visualization-only `$rendered_ci` field when `scaling` is `"matched"`
  or `"constant"`.

- col_participant, col_stimulus, col_response:

  Column names.

- method:

  Brief-RC variant: `"briefrc12"` (12 alternatives per trial, 6
  original + 6 inverted; the default) or `"briefrc20"` (20 alternatives
  per trial, 10 original + 10 inverted). Both variants are validated in
  Schmitz, Rougier, & Yzerbyt (2024). The CI computation is identical
  for both; the argument is kept as metadata so that downstream code and
  reports can record which paradigm produced the data.

- scaling:

  Visualization-only scaling for the optional `$rendered_ci` field. One
  of `"none"` (default), `"matched"` (stretch mask to base range, then
  add to base) or `"constant"` (multiply mask by `scaling_constant`,
  then add to base). The mathematical `$signal_matrix` is always the raw
  unscaled mask.

- scaling_constant:

  Numeric multiplier used when `scaling = "constant"`. Ignored
  otherwise.

- base_image_path:

  **Deprecated.** Use `base_image` (which accepts both a numeric matrix
  and a path). The old name still works for one release with a
  deprecation warning.

## Value

A list with `signal_matrix`, optionally `rendered_ci`, `participants`,
`img_dims`, `scaling`, and `method` (the Brief-RC variant the call was
made with).

## Details

Formula (Schmitz's `genMask()` exactly):

    X    <- data.table(response, stim)
    X    <- X[, .(response = mean(response)), stim]  # collapse duplicates
    mask <- (noiseMatrix[, X$stim] %*% X$response) / length(X$response)

The `length(X$response)` denominator is the number of unique pool ids
chosen by that participant, not the raw trial count. If a participant
chooses the same stimulus on two trials with opposite responses, those
two cancel.

The formula is symmetric in the per-trial split (6/6 for Brief-RC 12,
10/10 for Brief-RC 20), so the same code path handles both variants. The
`method` argument is recorded as provenance metadata and validated; it
does not branch the computation.

## Reading the result

- `$signal_matrix` is the raw mask per producer; pass this and only this
  to reliability metrics or any external infoVal computation.

- `$rendered_ci`, when present, is `base + scaling(mask)` per producer.
  **Visualization only.**

- `$participants` and `$img_dims` are convenience metadata.

## Common mistakes

- Passing the "expanded" multi-row-per-trial response format. Brief-RC
  data is **one row per trial**; `stimulus` is the chosen pool id and
  `response` is `+1` (original) or `-1` (inverted). See Schmitz et
  al. (2024) sec. 3.1.2.

- Using `$rendered_ci` for downstream stats; it exists only because the
  user often wants to save a PNG.

## References

Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the brief
reverse correlation: an improved tool to assess visual representations.
*European Journal of Social Psychology*.
[doi:10.1002/ejsp.3100](https://doi.org/10.1002/ejsp.3100)

## See also

[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md),
[`run_reliability()`](https://olivethree.github.io/rcisignal/reference/run_reliability.md),
[`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md),
[`group_ci()`](https://olivethree.github.io/rcisignal/reference/group_ci.md)
(stage 2, optional)

## Examples

``` r
if (FALSE) { # \dontrun{
# In your study, you replace the simulated inputs with your own:
#   - `sim$data` -> your responses data frame (CSV via read_responses())
#   - `sim$noise_matrix` -> a numeric matrix `n_pixels x pool_size`
#     (load via read_noise_matrix() from the OSF txt that came with
#     your Brief-RC stimuli).
# The `simulate_briefrc_data()` helper just produces the same shapes so
# the example below is runnable end-to-end.
sim <- simulate_briefrc_data(n_per_condition = 10, n_trials = 60, seed = 1)
res <- ci_from_responses_briefrc(
  sim$data,
  noise_matrix = sim$noise_matrix,
  base_image   = sim$base_face
)
dim(res$signal_matrix)   # n_pixels x n_producers
rel_split_half(res$signal_matrix, n_permutations = 200L, seed = 1)
} # }
```
