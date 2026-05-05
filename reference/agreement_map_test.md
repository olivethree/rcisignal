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

  Optional logical vector of length `nrow(signal_matrix)`; the test is
  computed on the masked pixel subset. Pixels outside the mask are
  returned as `NA_real_` per-pixel and `FALSE` in the significant mask.

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

- `$alpha`, `$n_permutations`, `$n_participants`, `$mask`.

## See also

[`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md),
[`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md),
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# In a real pipeline, signal_matrix comes from earlier steps:
#   signal_matrix <- ci_from_responses_briefrc(...)$signal_matrix
# For a self-contained demo we fabricate a small synthetic input:
n_pix  <- 32L * 32L
n_prod <- 20L
set.seed(1)
signal_matrix <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)

res <- agreement_map_test(signal_matrix,
                          n_permutations = 500L, seed = 1)
print(res)
} # }
```
