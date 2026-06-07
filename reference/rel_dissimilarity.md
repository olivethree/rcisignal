# Between-condition dissimilarity with bootstrap confidence intervals

Quantifies overall dissimilarity between two conditions' group-level
classification images. The primary statistic is Euclidean distance
between the two group-mean CIs, reported both raw and normalised by
`sqrt(n_pixels)` for cross-resolution comparability. Percentile
bootstrap 95% confidence intervals are computed by resampling
participants with replacement within each condition.

Pair with
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)
when you also want to know **where** the two conditions differ.

## Usage

``` r
rel_dissimilarity(
  signal_matrix_a,
  signal_matrix_b,
  paired = FALSE,
  n_boot = 2000L,
  ci_level = 0.95,
  null = c("none", "permutation"),
  n_permutations = 2000L,
  mask = NULL,
  seed = NULL,
  progress = TRUE,
  acknowledge_scaling = FALSE
)
```

## Arguments

- signal_matrix_a, signal_matrix_b:

  Pixels x participants, base-subtracted. Row counts must match.

- paired:

  Logical. `FALSE` (default) for between-subjects: participants
  resampled within A and B independently. `TRUE` for within-subjects: A
  and B share a single resample index per replicate so the paired
  covariance structure is preserved.

- n_boot:

  Bootstrap replicates. Default 2000.

- ci_level:

  Confidence level. Default 0.95.

- null:

  One of `"none"` (default) or `"permutation"`. Set `"permutation"` for
  any above-chance claim about the distance: it builds an empirical
  chance baseline by shuffling condition labels (or sign-flipping paired
  differences) and recomputing the distance. The default bootstrap CI is
  a precision interval only, not a test against zero (see the section on
  why a bootstrap CI is not a test against zero).

- n_permutations:

  Integer. Number of null iterations when `null = "permutation"`.
  Default 2000.

- mask:

  Optional logical vector of length `nrow(signal_matrix_a)`
  (column-major) restricting the Euclidean / correlation computation to
  a region. Build with
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  (parametric oval and sub-regions) or
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
  (PNG/JPEG mask).

- seed:

  Optional integer; RNG state restored on exit.

- progress:

  Show a `cli` progress bar.

- acknowledge_scaling:

  Logical. When `FALSE` (default), the shared `assert_raw_signal()`
  helper errors on a known-rendered matrix on either side.

## Value

Object of class `rcisignal_rel_dissim`.

## Details

For the observed statistics and each bootstrap replicate `i`:

    mean_a      = rowMeans(signal_matrix_a)
    mean_b      = rowMeans(signal_matrix_b)
    observed_dist            = sqrt(sum((mean_a - mean_b)^2))
    observed_dist_normalised = observed_dist / sqrt(n_pixels)

Percentile CI via base R
[`quantile()`](https://rdrr.io/r/stats/quantile.html); no `boot`
dependency.

## Why Euclidean and not Pearson correlation as the primary

Two base-subtracted CIs share systematic image-domain spatial structure
(face shape, oval signal support, low-frequency Gaussian-noise
smoothness) that pushes their correlation above zero even when the
underlying mental representations are unrelated. Absolute correlation
values therefore do not cleanly mean "these conditions are similar": two
unrelated traits can easily produce `r` well above zero from the shared
image scaffolding alone, so a "high" correlation does not by itself
license a similarity claim. Euclidean distance does not share this
baseline issue, which is why it is the primary statistic here.

## Why a bootstrap CI on the L2 distance is not a test against zero

The Euclidean distance is a non-negative L2 norm of the difference
between two group means. Resampling producers with replacement adds
variance to each resampled group mean, and that variance enters the
squared-distance sum at every pixel, so
`E[(mean_a_boot - mean_b_boot)^2]` is roughly
`(mean_a - mean_b)^2 + Var(mean_a_boot) + Var(mean_b_boot)`. Summed over
thousands of pixels the upward bias is large: the bootstrap distribution
of the distance sits *above* the observed value, and its confidence
interval almost always excludes zero even when the two conditions are
drawn from the same population. The bootstrap CI on this distance
therefore does not answer "is there a difference?". It answers only "how
stable is my distance estimate under producer resampling?", and it is
useful for comparing relative magnitudes across contrasts where the bias
is roughly constant.

To test whether the observed distance is larger than chance, set
`null = "permutation"`. The permutation null shuffles producer condition
labels (or sign-flips per-producer differences in paired designs),
recomputes the distance, and repeats. Its distribution is centred not at
zero but at the chance distance between any two random subgroups of
producers, which is itself positive because each subgroup mean carries
within-group noise. Read the observed distance against the upper tail of
*that* distribution via `$d_z`, `$d_ratio`, `$d_null_p95`, or a
permutation *p* computed as the proportion of `$null_distribution`
values at or above `$euclidean`. The Pearson fields carry the
mirror-image bias: resampling *attenuates* correlation, so the bootstrap
CI on `$correlation` sits *below* the observed value. Benchmark *r*
against a permutation null as well, never against zero.

The Pearson correlation fields (`$correlation`, `$boot_cor`, `$ci_cor`,
`$boot_se_cor`) are returned alongside it as a secondary summary, for
users who need to report `r` for comparability with prior literature or
with another analysis pipeline. They are not recommended as a standalone
similarity score. If you do report `r`, interpret it **relatively** (the
ordering of `r` across multiple condition pairs is more defensible than
any single absolute value) and **benchmark against a permutation null**
rather than against zero. The image-domain scaffolding shifts the chance
baseline upward; "above zero" is not the right comparison.

## Reading the plot

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the returned
object renders the bootstrap distributions as two side-by-side
histograms (Euclidean distance on the left, Pearson r on the right). The
Pearson panel is rendered in grey because it is a secondary summary, not
recommended as a standalone similarity score; see the description for
why and how to use it carefully if needed.

- The shaded vertical band marks the percentile CI at `ci_level`.

- The vertical line marks the *observed* statistic on the real data (not
  a bootstrap mean).

- The Euclidean CI band is a *precision* interval on the distance
  estimate, not a test against zero. Producer resampling inflates the
  distance, so this band sits above the observed value and almost always
  excludes zero even for identical conditions; do not read "excludes
  zero" as evidence of a difference. For an above-chance test, run with
  `null = "permutation"`. The numbers are returned in `$ci_dist`.

- For visual comparison across multiple contrasts, pass each
  `rel_dissimilarity()` result to
  [`plot_dissimilarity_grid()`](https://olivethree.github.io/rcisignal/reference/plot_dissimilarity_grid.md),
  which lays them out as labelled CI bars on a shared axis.

- For a *spatial* picture of where the two conditions differ, pair this
  with
  [`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)
  (or use
  [`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md)
  to run both in one call).

## Reading the result

- `$euclidean`, observed Euclidean distance between group means (primary
  statistic).

- `$euclidean_normalised`, `$euclidean / sqrt(n_pixels)`. Use for
  cross-resolution comparisons.

- `$boot_dist`, `$ci_dist`, `$boot_se_dist`, bootstrap distribution,
  percentile CI, and SE of the Euclidean distance. The CI is a precision
  interval on the estimate, not a test against zero: it is biased upward
  by resampling and excludes zero even for identical conditions. Use
  `null = "permutation"` to test against chance.

- `$null` (character), the null mode used.

- `$null_distribution`, when `null != "none"`: numeric vector of
  per-iteration Euclidean distances under the chosen null.

- `$d_null_p95`, 95th percentile of the null distribution.

- `$d_z`, z-equivalent effect size:
  `(observed_d - mean(null)) / sd(null)`.

- `$d_ratio`, observed Euclidean over the null median.

- `$correlation`, `$boot_cor`, `$ci_cor`, `$boot_se_cor`: Pearson
  correlation of the group means and its bootstrap summaries. Secondary;
  not recommended as a standalone similarity score (see the "Why
  Euclidean" section). If reporting, prefer relative comparisons across
  pairs against a permutation null.

- `$n_boot`, `$ci_level`, `$paired`, metadata.

## References

Efron, B., & Tibshirani, R. J. (1994). *An introduction to the
bootstrap*. Chapman & Hall / CRC.

## See also

[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md),
[`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal call-signature demo with two synthetic inputs.
n_pix  <- 32L * 32L
n_prod <- 20L
set.seed(1)
signal_matrix_a <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
signal_matrix_b <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
rel_dissimilarity(signal_matrix_a, signal_matrix_b,
                  n_boot = 200L, seed = 1)
} # }

if (FALSE) { # \dontrun{
# Same function, richer input: signal planted in different face regions
# (eyes vs mouth). The Euclidean distance and its bootstrap CI should
# be well above zero, reflecting genuine spatial divergence.
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
d <- rel_dissimilarity(sig_eyes, sig_mouth, n_boot = 500L, seed = 1)
# Bootstrap distribution + observed Euclidean + 95% CI band.
plot(d, main = "Eyes vs Mouth: bootstrap dissimilarity")
} # }
```
