# Leave-one-out influence screening

Screens for influential producers by asking: "if this producer is
removed from the sample, how much does the group classification image
change?" Producers whose removal moves the group pattern
disproportionately are flagged for inspection. This is an **influence /
outlier screening** tool, not a reliability statistic.

## Usage

``` r
rel_loo(
  signal_matrix,
  flag_threshold = 2.5,
  flag_method = c("mad", "sd"),
  flag_threshold_sd = NULL,
  mask = NULL
)
```

## Arguments

- signal_matrix:

  Pixels x participants, base-subtracted.

- flag_threshold:

  Numeric multiplier on `sd` (or `mad`) below the centre. Default 2.5.

- flag_method:

  One of `"mad"` (default) or `"sd"`. SD/mean is retained for backwards
  compatibility and emits a one-time per-session deprecation message.

- flag_threshold_sd:

  Deprecated alias for `flag_threshold`.

- mask:

  Optional logical vector of length `nrow(signal_matrix)` (column-major)
  restricting the LOO correlation to a region. Build with
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  (parametric oval and sub-regions) or
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
  (PNG/JPEG mask).

## Value

Object of class `rcisignal_rel_loo`.

## Details

For each producer `i`, compute the Pearson correlation between the
full-sample group CI and the group CI recomputed without that producer:

    full        <- rowMeans(signal_matrix)
    r_loo[i]    <- cor(full, rowMeans(signal_matrix[, -i]))

Because the full-sample mean and the leave-one-out mean share
`(N - 1) / N` of their data, `r_loo` values are near 1 by construction
even on noisy data, typically `[0.95, 0.999]` at `N = 30`. The absolute
level of `r_loo` is not informative. What is informative is the relative
ordering: producers whose `r_loo` sits clearly below the pack are
candidates for inspection. The function therefore returns a z-scored
version of `r_loo` in `$z_scores` using the same centre / spread
estimators as the flagging rule. `$z_scores` is the recommended quantity
to plot or report.

Two flagging rules:

- `"mad"` (default): flag producers with
  `r_loo < median(r) - flag_threshold * mad(r)`. Robust to the few
  atypical producers RC datasets often contain.

- `"sd"`: flag producers with
  `r_loo < mean(r) - flag_threshold * sd(r)`. Sensitive to the very
  outliers it is trying to detect; emits a one-time per-session
  deprecation message and is scheduled for removal in v0.2.0.

Default `flag_threshold = 2.5` is calibrated so a 30-producer dataset
flags roughly 0.3 producers by chance under `"sd"`, rather than the ~1.5
a 2-SD rule would produce. Under `"mad"` this is roughly comparable
thanks to MAD's 1.4826 consistency factor.

## What this function is, and is not

`rel_loo()` is an influence-screening diagnostic. It answers "which
producers disproportionately shape the group CI?", not "how reliable is
the group CI?". For reliability use
[`rel_split_half()`](https://olivethree.github.io/rcisignal/reference/rel_split_half.md)
or
[`rel_icc()`](https://olivethree.github.io/rcisignal/reference/rel_icc.md).
A flag does not mean the producer is "bad"; it means the producer's
individual CI sits far enough from the group pattern that the data
deserve a second look.

## Reading the result

- `$z_scores`, named numeric vector, per-producer standardized
  influence. The recommended quantity to plot or threshold.

- `$correlations`, named numeric vector, raw per-producer `r_loo`
  values.

- `$mean_r`, `$sd_r`, `$median_r`, `$mad_r`, centre / spread under each
  rule.

- `$threshold`, raw cutoff value on `r_loo`.

- `$flagged`, character vector of producer ids below threshold.

- `$summary_df`, one row per producer with `correlation`, `z_score`, and
  `flag`, sorted by `z_score`.

- `$flag_method`, `$flag_threshold`.

## Common mistakes

- Reading `r_loo` as a reliability. An `r_loo` of .98 does not mean the
  CI is 98% reliable; it means a single producer's removal changed the
  group mean by 2%.

- Treating `$flagged` as "drop these producers". Investigate first;
  cross-check with the rcisignal input-side
  [`run_diagnostics()`](https://olivethree.github.io/rcisignal/reference/run_diagnostics.md)
  to rule out response-coding errors.

- Lowering `flag_threshold` below 2 to flag more producers; that trades
  real signal for noise. Use `flag_method = "mad"` instead if the SD
  rule is dominated by outliers.

## Reliability metrics expect raw masks

Operates on the raw mask; results may be distorted if `signal_matrix`
was extracted from rendered (scaled) PNGs.

## See also

[`rel_loo_z()`](https://olivethree.github.io/rcisignal/reference/rel_loo_z.md)
for a tidy z-score accessor;
[`rel_split_half()`](https://olivethree.github.io/rcisignal/reference/rel_split_half.md),
[`rel_icc()`](https://olivethree.github.io/rcisignal/reference/rel_icc.md)
for reliability metrics proper;
[`run_reliability()`](https://olivethree.github.io/rcisignal/reference/run_reliability.md).

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

r <- rel_loo(signal_matrix)
print(r)
} # }
```
