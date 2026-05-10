# Per-producer informational value with trial-count-matched reference

Computes a z-scored informational value (infoVal) for each producer's
classification image, using a reference distribution matched to that
producer's trial count. Handles both 2IFC and Brief-RC paradigms with a
single function: the difference is entirely in what the user passes as
`noise_matrix`, not in how the statistic is computed.

## Usage

``` r
infoval(
  signal_matrix,
  noise_matrix,
  trial_counts,
  iter = 10000L,
  mask = NULL,
  with_replacement = c("auto", TRUE, FALSE),
  cache_path = NULL,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- signal_matrix:

  Pixels x participants numeric matrix of raw masks (as returned by
  [`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md)
  or
  [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md)).

- noise_matrix:

  Pixels x pool-size numeric matrix of noise patterns. Row count must
  match `signal_matrix`.

- trial_counts:

  Named integer vector of trial counts per producer. Names must match
  `colnames(signal_matrix)`.

- iter:

  Reference-distribution Monte Carlo size. Default 10000.

- mask:

  Optional logical vector of length `nrow(signal_matrix)`
  (column-major). When supplied, both observed and reference norms are
  computed on the masked pixel subset. Build with
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  (parametric oval and sub-regions) or
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
  (PNG/JPEG mask).

- with_replacement:

  Sampling regime for the reference distribution at the across-trials
  level (how stimulus ids are drawn when simulating a random producer).
  One of `"auto"` (default), `TRUE`, or `FALSE`. `"auto"` matches the
  standard convention: without replacement when a producer's trial count
  fits in the pool, with replacement otherwise. Set explicitly only when
  your design departs from this convention.

- cache_path:

  Optional path to an `.rds` file. When set and the file exists,
  reference distributions are loaded if the cache's `iter`, `n_pool`,
  and mask signature match; otherwise the reference is recomputed and
  saved to `cache_path`. The cache holds reference norms only, never
  observed data.

- seed:

  Optional integer; RNG state restored on exit.

- progress:

  Show a `cli` progress bar.

## Value

Object of class `rcisignal_rel_infoval` with fields described in
**Reading the result** above.

## Details

The observed statistic is the Frobenius norm of producer j's
classification image mask, optionally restricted to a logical `mask`
(e.g. an oval face region via
[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)):

    norm_j = sqrt(sum( signal_matrix[mask, j]^2 ))

The reference distribution for producer j is built by simulating `iter`
random masks at the same trial count `trial_counts[j]`, mirroring
Schmitz, Rougier, & Yzerbyt (2024)'s `genMask()` construction:

1.  Sample `trial_counts[j]` stimulus ids uniformly from `1:n_pool`
    (without replacement when `trial_counts[j] <= n_pool`).

2.  Sample `trial_counts[j]` responses uniformly from `{-1, +1}` (always
    with replacement).

3.  Collapse `response` by `stim` via
    [`mean()`](https://rdrr.io/r/base/mean.html).

4.  Build the mask:
    `noise_matrix[, unique_stims] %*% mean_response / n_unique_stims`.

5.  Apply `mask` (if supplied) and compute the Frobenius norm.

Per-producer z-score: `(norm_j - median(ref)) / mad(ref)`. Producers
sharing a trial count share a reference; this makes the simulation
efficient.

Choice of statistic: the **Frobenius norm** is the correct form per
Schmitz, Rougier, & Yzerbyt (2020 comment + 2020 erratum). The original
[`rcicr::computeInfoVal2IFC()`](https://rdrr.io/pkg/rcicr/man/computeInfoVal2IFC.html)
already uses it (`norm(matrix(target_ci[["ci"]]), "f")`); this function
matches that convention for both 2IFC and Brief-RC data.

Trial-count matching: the original
[`rcicr::generateReferenceDistribution2IFC()`](https://rdrr.io/pkg/rcicr/man/generateReferenceDistribution2IFC.html)
builds the reference using the full pool size. For Brief-RC the
producer's recorded `n_trials` is typically smaller than `n_pool`, so a
pool-size reference biases infoVal downward. `infoval()` closes this gap
by keying the reference on the actual per-producer recorded trial count.

What infoVal measures (and what it does not). The Frobenius norm is a
**magnitude** statistic, not a target-alignment one. Cross-paradigm
comparisons need care; stability and discriminability are different
questions and call for
[`rel_split_half()`](https://olivethree.github.io/rcisignal/reference/rel_split_half.md)
/
[`rel_icc()`](https://olivethree.github.io/rcisignal/reference/rel_icc.md)
/
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)
/
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md).

## References

Brinkman, L., Goffin, S., van de Schoot, R., van Haren, N. E. M.,
Dotsch, R., & Aarts, H. (2019). Quantifying the informational value of
classification images. *Behavior Research Methods*, 51(5), 2059-2073.
[doi:10.3758/s13428-019-01232-2](https://doi.org/10.3758/s13428-019-01232-2)

Schmitz, M., Rougier, M., & Yzerbyt, V. (2020). Comment on "Quantifying
the informational value of classification images": A miscomputation of
the infoVal metric. *Behavior Research Methods*, 52(3), 1383-1386.
[doi:10.3758/s13428-019-01295-1](https://doi.org/10.3758/s13428-019-01295-1)

Schmitz, M., Rougier, M., Yzerbyt, V., Brinkman, L., & Dotsch, R.
(2020). Erratum to: Comment on "Quantifying the informational value of
classification images". *Behavior Research Methods*, 52(4), 1800-1801.
[doi:10.3758/s13428-020-01367-7](https://doi.org/10.3758/s13428-020-01367-7)

Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the brief
reverse correlation: an improved tool to assess visual representations.
*European Journal of Social Psychology*.
[doi:10.1002/ejsp.3100](https://doi.org/10.1002/ejsp.3100)

## See also

[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md),
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md),
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal call-signature demo with three synthetic inputs.
n_pix  <- 32L * 32L
n_prod <- 20L
set.seed(1)
signal_matrix <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
noise_matrix  <- matrix(rnorm(n_pix * 4096L),  n_pix, 4096L)
trial_counts  <- rep(300L, n_prod)
iv <- infoval(signal_matrix, noise_matrix, trial_counts,
              iter = 200L, seed = 1)
print(iv)
} # }

if (FALSE) { # \dontrun{
# Same function, richer input: plant a strong signal in the eye region
# so per-producer infoVal z values mostly clear the 1.96 threshold.
sim <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "target",
  signal_region = "eyes", signal_strength = "strong", seed = 1
)
cis <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
trial_counts <- as.integer(table(sim$data$participant_id))
infoval(cis$signal_matrix, sim$noise_matrix, trial_counts,
        iter = 500L, seed = 1)
} # }
```
