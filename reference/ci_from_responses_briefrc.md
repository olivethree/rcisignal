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
  base_image_path,
  participant_col = "participant_id",
  stimulus_col = "stimulus",
  response_col = "response",
  method = c("briefrc12", "briefrc20"),
  scaling = c("none", "matched", "constant"),
  scaling_constant = NULL
)
```

## Arguments

- responses:

  Data frame with one row per trial. Must contain the columns named by
  `participant_col`, `stimulus_col`, `response_col`. `response` values
  must be in `{-1, +1}`.

- rdata_path, noise_matrix:

  Exactly one must be supplied. Provide `rdata_path` to read the noise
  matrix from an rcicr `.Rdata`, or pass a pre-loaded `noise_matrix`
  directly.

- base_image_path:

  Path to the base face image (PNG / JPEG). Used to validate the noise
  matrix shape and (when `scaling` is not `"none"`) to render the
  visualisation-only `$rendered_ci` field.

- participant_col, stimulus_col, response_col:

  Column names.

- method:

  Brief-RC variant: `"briefrc12"` (12 alternatives per trial, 6
  original + 6 inverted; the default) or `"briefrc20"` (20 alternatives
  per trial, 10 original + 10 inverted). Both variants are validated in
  Schmitz, Rougier, & Yzerbyt (2024). The CI computation is identical
  for both; the argument is kept as metadata so that downstream code and
  reports can record which paradigm produced the data.

- scaling:

  Visualisation-only scaling for the optional `$rendered_ci` field. One
  of `"none"` (default), `"matched"` (stretch mask to base range, then
  add to base) or `"constant"` (multiply mask by `scaling_constant`,
  then add to base). The mathematical `$signal_matrix` is always the raw
  unscaled mask.

- scaling_constant:

  Numeric multiplier used when `scaling = "constant"`. Ignored
  otherwise.

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
  **Visualisation only.**

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
[`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- ci_from_responses_briefrc(
  responses       = my_responses,
  rdata_path      = "data/rcicr_stimuli.Rdata",
  base_image_path = "data/base.jpg"
)
rel_split_half(res$signal_matrix, seed = 1L)
} # }
```
