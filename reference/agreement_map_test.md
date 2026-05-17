# Per-pixel inferential agreement map with FWER control

Within a single condition, tests at each pixel whether the
producer-level signal differs from zero (one-sample t). The permutation
null is built by random sign-flip per producer (exact under the
assumption that, under H0, the producer's signal contribution is
symmetric around zero). Family-wise error is controlled across pixels by
the maximum \|t\| statistic.

Use this when you want a per-pixel inferential overlay on a descriptive
agreement-map plot, typically paired with
[`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)
so the significance contours are rendered on top of the observed group
CI.

## Usage

``` r
agreement_map_test(
  signal_matrix,
  n_permutations = 5000L,
  alpha = 0.05,
  mask = NULL,
  seed = NULL,
  progress = TRUE,
  acknowledge_scaling = FALSE
)
```

## Arguments

- signal_matrix:

  Pixels x participants, base-subtracted.

- n_permutations:

  Integer. Number of sign-flip iterations. Default 5000.

- alpha:

  Numeric. Significance level. Default 0.05.

- mask:

  Optional logical vector of length `nrow(signal_matrix)`
  (column-major). When supplied, the test is computed on the masked
  pixel subset. Pixels outside the mask are returned as `NA_real_`
  per-pixel and `FALSE` in the significant mask. Build with
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  (parametric oval and sub-regions) or
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
  (PNG/JPEG mask).

- seed:

  Optional integer.

- progress:

  Show a `cli` progress bar.

- acknowledge_scaling:

  Logical. Forwarded to `assert_raw_signal()`.

## Value

Object of class `rcisignal_rel_agreement_map_test` with:

- `$observed_t`: per-pixel one-sample t.

- `$pmap`: per-pixel p-value under the max-\|t\| null.

- `$significant_mask`: logical, `pmap < alpha`.

- `$null_distribution`: numeric vector of `max_abs_t` per permutation.

- `$img_dims`: `c(nrow, ncol)` inferred from
  `attr(signal_matrix, "img_dims")` or, when absent, from a square pixel
  count. Used by the S3
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method to
  reshape the per-pixel vectors back to an image.

- `$alpha`, `$n_permutations`, `$n_participants`, `$mask`.

## Reading the result and the paired plot

- `$observed_t`: per-pixel one-sample t (sign indicates the direction of
  producer agreement; magnitude indicates strength).

- `$pmap`: per-pixel p under the max-\|t\| null. Already FWE-corrected
  at the alpha you chose; no further adjustment needed.

- `$significant_mask`: logical vector of pixels with `pmap < alpha`.
  This is the field that `plot_ci_overlay(..., test = result)` traces
  with black contours on top of the group CI overlay.

- Pass the result to
  `plot_ci_overlay(signal_matrix, base_image, test = result)` to get the
  canonical figure: blue = positive producer-mean signal, red =
  negative, opacity = magnitude, black contours = FWE-significant pixels
  at this test's alpha. To plot only the inferential mask (no continuous
  overlay), pass
  `signal_matrix = result$observed_t * result$significant_mask` or use
  `mask = result$significant_mask` to clip the overlay to the
  significant region.

## See also

[`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md),
[`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md),
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal call-signature demo with a synthetic input.
n_pix  <- 32L * 32L
n_prod <- 20L
set.seed(1)
signal_matrix <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
res <- agreement_map_test(signal_matrix,
                          n_permutations = 500L, seed = 1)
print(res)
} # }

if (FALSE) { # \dontrun{
# Same function, richer input: plant a strong signal in the eye region
# so the FWER-controlled test has something to detect.
sim <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "target",
  signal_region = "eyes", signal_strength = "strong", seed = 1
)
cis <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
agreement_map_test(cis$signal_matrix, n_permutations = 500L, seed = 1)
} # }

if (FALSE) { # \dontrun{
# Canonical pairing: feed the test result to plot_ci_overlay() so the
# FWE-significant pixels are outlined in black on top of the group CI.
sim   <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "target",
  signal_region = "eyes", signal_strength = "strong", seed = 1
)
cis   <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
agree <- agreement_map_test(cis$signal_matrix,
                            n_permutations = 500L, seed = 1)
plot_ci_overlay(cis$signal_matrix, sim$base_face,
                test = agree,
                main = "CI overlay + FWE-significant pixels")
} # }
```
