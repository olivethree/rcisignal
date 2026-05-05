# Intraclass correlation coefficients via direct mean squares

Reports the reliability of producer-level CI signals as a single scalar
(or two scalars) per condition. Use **ICC(3,1)** to ask "how informative
is one producer's CI as a noisy estimate of the group pattern?" and
**ICC(3,k)** to ask "how stable is the group-mean CI this experiment
produced?". Most papers report ICC(3,k) as the headline.

## Usage

``` r
rel_icc(
  signal_matrix,
  variants = c("3_1", "3_k"),
  mask = NULL,
  acknowledge_scaling = FALSE
)
```

## Arguments

- signal_matrix:

  Pixels x participants (targets x raters), base-subtracted.

- variants:

  Character vector of which ICC variants to return. Subset of
  `c("3_1", "3_k", "2_1", "2_k")`. Defaults to `c("3_1", "3_k")`.

- mask:

  Optional logical vector of length `nrow(signal_matrix)` restricting
  computation to a region (e.g., from
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  or
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)).

- acknowledge_scaling:

  Logical. When `FALSE` (default), the shared `assert_raw_signal()`
  helper errors on a known-rendered matrix.

## Value

Object of class `rcisignal_rel_icc`.

## Details

Reports **ICC(3,1)** and **ICC(3,k)**, two-way mixed model with pixels
fixed and participants random. Pixels are a fixed `img_size x img_size`
grid (not a random sample), so ICC(2,\*) is mis-specified even when
numerically similar. ICC(2,1) and ICC(2,k) are available via `variants`
for comparability with reports that use the two-way-random model.

Computed directly from ANOVA mean squares (never via
[`psych::ICC()`](https://rdrr.io/pkg/psych/man/ICC.html), which
allocates intermediates that blow memory on a 262,144 x 30 matrix).

## What this ICC is, and is not

Most ICCs reported in the reverse-correlation literature are
trait-rating reliability (phase-2 raters scoring CIs on trait
dimensions). rcisignal's ICC is structurally different: it operates on
the pixel-level signal produced by the original producers. No phase-2
rating study is involved. This sidesteps the two-phase design Cone,
Brown-Iannuzzi, Lei, & Dotsch (2021) showed inflates Type I error.

## Reading the result

- `$icc_3_1`, `$icc_3_k`, two-way-mixed ICCs for single rater and
  average rater respectively. `$icc_2_1` / `$icc_2_k` only present if
  requested via `variants`.

- `$ms_rows`, `$ms_cols`, `$ms_error`, the underlying ANOVA mean
  squares.

- `$n_raters` (= participants), `$n_targets` (= pixels), `$model`,
  `$variants`.

## Common mistakes

- Reporting ICC(2,*) as the headline. ICC(2,*) assumes targets are a
  random sample from a target population (the right model when raters
  score a random subset of stimuli drawn from a larger pool). Here the
  "targets" are the same fixed `img_size x img_size` pixel grid in every
  CI, not a sample, so the random-targets assumption does not hold.
  Numerically close to ICC(3,\*) at high pixel counts because the bias
  from mis-specification shrinks with target count, but the model is
  wrong for this object.

- Comparing this ICC to a phase-2 trait-rating ICC from a different
  paper. Different statistical objects.

## Reliability metrics expect raw masks

ICC is variance-based and strongly sensitive to any scaling step. Inputs
with `attr(., "source") == "rendered"` (set automatically by Mode 1
readers like
[`extract_signal()`](https://olivethree.github.io/rcisignal/reference/extract_signal.md))
error unless `acknowledge_scaling = TRUE`.

## References

Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: uses in
assessing rater reliability. *Psychological Bulletin*, 86(2), 420-428.
[doi:10.1037/0033-2909.86.2.420](https://doi.org/10.1037/0033-2909.86.2.420)

McGraw, K. O., & Wong, S. P. (1996). Forming inferences about some
intraclass correlation coefficients. *Psychological Methods*, 1(1),
30-46.
[doi:10.1037/1082-989X.1.1.30](https://doi.org/10.1037/1082-989X.1.1.30)

Cone, J., Brown-Iannuzzi, J. L., Lei, R., & Dotsch, R. (2021). Type I
error is inflated in the two-phase reverse correlation procedure.
*Social Psychological and Personality Science*, 12(5), 760-768.
[doi:10.1177/1948550620938616](https://doi.org/10.1177/1948550620938616)

## See also

[`rel_split_half()`](https://olivethree.github.io/rcisignal/reference/rel_split_half.md),
[`rel_loo()`](https://olivethree.github.io/rcisignal/reference/rel_loo.md),
[`run_reliability()`](https://olivethree.github.io/rcisignal/reference/run_reliability.md)

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

r <- rel_icc(signal_matrix)
print(r)
} # }
```
