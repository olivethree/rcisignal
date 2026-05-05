# Permuted split-half reliability with Spearman-Brown correction

Measures how stable a condition's group-level CI is across random halves
of its producer set. For each of `n_permutations` iterations: randomly
partition the producers (columns of `signal_matrix`) into two
roughly-equal halves, average each half to get two group vectors, and
compute the Pearson correlation.

Use this to answer "would I get the same group CI if I'd run a different
half of my participants?".

For odd N, one randomly-chosen producer is dropped per permutation
(re-drawn each iteration, not fixed) so both halves have `floor(N/2)`
producers.

Reports both the mean per-permutation `r` and the Spearman-Brown
corrected reliability `r_SB = (2 r) / (1 + r)`. 95% CI is the 2.5 / 97.5
percentile of the permutation distribution.

## Usage

``` r
rel_split_half(
  signal_matrix,
  n_permutations = 2000L,
  null = c("none", "permutation", "random_responders"),
  noise_matrix = NULL,
  mask = NULL,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- signal_matrix:

  Pixels x participants, base-subtracted.

- n_permutations:

  Integer number of random splits. Default 2000.

- null:

  One of `"none"` (default), `"permutation"`, or `"random_responders"`.

  - `"none"` skips the empirical-null computation; `$null_distribution`
    is `NULL`.

  - `"permutation"`: at each iteration, generate fresh Gaussian noise
    per producer (no shared spatial structure); recompute `r_hh` per
    iteration. The empirical distribution becomes the chance baseline.

  - `"random_responders"`: simulate `ncol(signal_matrix)` producers
    responding at chance using `genMask()` machinery on the supplied
    `noise_matrix`. Closer to the empirical chance baseline of an actual
    RC experiment.

- noise_matrix:

  Pixels x pool-size numeric matrix of noise patterns. Required when
  `null = "random_responders"`.

- mask:

  Optional logical vector of length `nrow(signal_matrix)` restricting
  computation to a region of the image (e.g., from
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  or
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)).

- seed:

  Optional integer; if set, results are reproducible. The caller's
  global RNG state is restored on exit.

- progress:

  Show a `cli` progress bar.

## Value

An object of class `rcisignal_rel_split_half`. Fields described in
**Reading the result** above.

## Reading the result

- `$r_hh`, mean per-permutation Pearson r between the two halves.

- `$r_sb`, Spearman-Brown projected reliability of the full sample. This
  is usually the headline number to report.

- `$ci_95`, `$ci_95_sb`, percentile 95% CIs on each.

- `$distribution`, the full per-permutation r vector for plotting or
  custom CIs.

- `$null` (character), the null mode used.

- `$null_distribution`, when `null != "none"`: numeric vector of
  per-iteration `r_hh` under the chosen null. `NULL` otherwise.

- `$r_hh_null_p95`, 95th percentile of the null distribution. `NA_real_`
  when `null = "none"`.

- `$r_hh_excess`, observed `r_hh` minus the null median.

- `$r_sb_excess`, the same excess applied to Spearman-Brown projected
  `r_sb`.

- `$n_participants`, `$n_permutations`, metadata.

## Common mistakes

- Reporting `$r_hh` instead of `$r_sb` when the audience cares about the
  full-sample CI's reliability (almost always).

- Running with a small `n_permutations` (\< 200) to save time and then
  trusting the CI bounds; CIs widen with low n_perm.

## Reliability metrics expect raw masks

Operates on the raw mask; results may be distorted if `signal_matrix`
was extracted from rendered (scaled) PNGs. See the rcisignal
output-reliability vignette chapter on raw vs. rendered.

## References

Brinkman, L., Todorov, A., & Dotsch, R. (2017). Visualising mental
representations: A primer on noise-based reverse correlation in social
psychology. *European Review of Social Psychology*, 28(1), 333-361.
[doi:10.1080/10463283.2017.1381469](https://doi.org/10.1080/10463283.2017.1381469)

Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: uses in
assessing rater reliability. *Psychological Bulletin*, 86(2), 420-428.
[doi:10.1037/0033-2909.86.2.420](https://doi.org/10.1037/0033-2909.86.2.420)

## See also

[`rel_loo()`](https://olivethree.github.io/rcisignal/reference/rel_loo.md),
[`rel_icc()`](https://olivethree.github.io/rcisignal/reference/rel_icc.md),
[`run_reliability()`](https://olivethree.github.io/rcisignal/reference/run_reliability.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal call-signature demo with a synthetic input. Split-half
# correlation will hover near zero because the input is pure noise.
n_pix  <- 32L * 32L
n_prod <- 20L
set.seed(1)
signal_matrix <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
rel_split_half(signal_matrix, n_permutations = 200L, seed = 1)
} # }

if (FALSE) { # \dontrun{
# Same function, richer input: plant a strong signal in the eye region.
# Split-half correlation should now be measurably positive.
sim <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "target",
  signal_region = "eyes", signal_strength = "strong", seed = 1
)
cis <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
rel_split_half(cis$signal_matrix, n_permutations = 200L, seed = 1)
} # }
```
