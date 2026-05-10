# Vectorised pixel-wise t-test (independent or paired)

At every pixel, tests whether condition A's mean signal differs from
condition B's. Two modes:

- `paired = FALSE` (default): independent-samples Welch t per pixel.
  Correct when producers in A and B are different people
  (between-subjects design).

- `paired = TRUE`: paired t per pixel on the per-producer difference
  `A - B`. Correct when the same producers contributed to both
  conditions (within-subjects design).

Returns a numeric vector of t-values, length `n_pixels`. Used as an
intermediate by
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md);
not intended as a standalone inferential test.

## Usage

``` r
pixel_t_test(
  signal_matrix_a,
  signal_matrix_b,
  paired = FALSE,
  mask = NULL,
  acknowledge_scaling = FALSE
)
```

## Arguments

- signal_matrix_a, signal_matrix_b:

  Pixels x participants, base-subtracted. Row counts must match. When
  `paired = TRUE` the column counts must also match, and column names
  must correspond to the same producer across matrices.

- paired:

  Logical. `FALSE` (default) uses independent Welch t; `TRUE` uses
  paired t.

- mask:

  Optional logical vector of length `nrow(signal_matrix_a)`
  (column-major). Both matrices are subsetted with the same mask before
  computing t; the returned vector is then of length `sum(mask)`. Build
  with
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  (parametric oval and sub-regions) or
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
  (PNG/JPEG mask).

- acknowledge_scaling:

  Logical. When `FALSE` (default), the shared `assert_raw_signal()`
  helper errors on a known-rendered matrix. Cascades to internal
  `pixel_t_test()` calls from
  [`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md).

## Value

Numeric vector of length `nrow(signal_matrix_a)` (or `sum(mask)` if
`mask` is supplied). Pixels with zero variance get `0` rather than
`NaN`.

## Reliability metrics expect raw masks

Welch t and paired t are variance-based and sensitive to scaling. Inputs
with `attr(., "source") == "rendered"` (set automatically by Mode 1
readers like
[`extract_signal()`](https://olivethree.github.io/rcisignal/reference/extract_signal.md))
error unless `acknowledge_scaling = TRUE`.

## See also

[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal call-signature demo with two synthetic inputs.
n_pix  <- 32L * 32L
n_prod <- 20L
set.seed(1)
signal_matrix_a <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
signal_matrix_b <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
summary(pixel_t_test(signal_matrix_a, signal_matrix_b))
} # }

if (FALSE) { # \dontrun{
# Same function, richer input: simulate two conditions with signal
# planted in different face regions (eyes vs mouth). The pixel-wise
# t-map should be large-positive around the eyes and large-negative
# around the mouth — pixels where the conditions disagree.
sim_eyes  <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "x",
  signal_region = "eyes", signal_strength = "strong", seed = 1
)
sim_mouth <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "x",
  signal_region = "mouth", signal_strength = "strong", seed = 2
)
sig_eyes  <- ci_from_responses_briefrc(
  sim_eyes$data, noise_matrix = sim_eyes$noise_matrix)$signal_matrix
sig_mouth <- ci_from_responses_briefrc(
  sim_mouth$data, noise_matrix = sim_mouth$noise_matrix)$signal_matrix
summary(pixel_t_test(sig_eyes, sig_mouth))
} # }
```
