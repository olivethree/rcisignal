# Cluster-based permutation test with family-wise error control

Pixel-level discriminability test between two condition signal matrices.
Returns either (a) the classical threshold-based cluster test (Maris &
Oostenveld, 2007) with cluster mass as the statistic, or (b)
threshold-free cluster enhancement (TFCE; Smith & Nichols, 2009), which
integrates across thresholds and removes the arbitrary cluster-forming
cutoff. Both paths use the same stratified permutation scheme and the
maximum-statistic approach to family-wise error control (Nichols &
Holmes, 2002).

Pair with
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
when you also want an overall magnitude summary.

## Usage

``` r
rel_cluster_test(
  signal_matrix_a,
  signal_matrix_b,
  img_dims = NULL,
  method = c("threshold", "tfce"),
  paired = FALSE,
  n_permutations = 2000L,
  cluster_threshold = 2,
  tfce_H = 2,
  tfce_E = 0.5,
  tfce_n_steps = 100L,
  alpha = 0.05,
  mask = NULL,
  seed = NULL,
  progress = TRUE,
  acknowledge_scaling = FALSE
)
```

## Arguments

- signal_matrix_a, signal_matrix_b:

  Pixels x participants, base-subtracted.

- img_dims:

  Integer `c(nrow, ncol)`. Inferred from
  `attr(signal_matrix_a, "img_dims")` when `NULL`.

- method:

  One of `"threshold"` (default) or `"tfce"`.

- paired:

  Logical. `FALSE` (default) for between-subjects: Welch t and
  stratified label permutation. `TRUE` for within-subjects: paired t on
  `A - B` and random sign-flip of matched pairs for the null.

- n_permutations:

  Default 2000.

- cluster_threshold:

  Absolute t-value threshold for forming clusters under
  `method = "threshold"`. Default 2.0.

- tfce_H, tfce_E:

  TFCE height and extent exponents. Defaults match Smith & Nichols
  (2009): 2.0 and 0.5.

- tfce_n_steps:

  Number of thresholds in the TFCE integration grid. Default 100.

- alpha:

  Significance level. Default 0.05.

- mask:

  Optional logical vector of length `n_pixels` (column-major)
  restricting cluster / TFCE inference to a region. The implementation
  uses a zero-out pattern (not row-subsetting) so the 2D image structure
  is preserved for 4-connectivity. Build with
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  (parametric oval and sub-regions) or
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
  (PNG/JPEG mask).

- seed:

  Optional integer; RNG state restored on exit.

- progress:

  Show a `cli` progress bar.

- acknowledge_scaling:

  Logical. When `FALSE` (default), the internal
  [`pixel_t_test()`](https://olivethree.github.io/rcisignal/reference/pixel_t_test.md)
  errors on a known-rendered matrix.

## Value

Object of class `rcisignal_rel_cluster_test`.

## Details

Common scaffold (both methods):

1.  Compute observed pixel-wise Welch t (see
    [`pixel_t_test()`](https://olivethree.github.io/rcisignal/reference/pixel_t_test.md)).

2.  Build the null via `n_permutations` stratified label permutations,
    each preserving `(N_A, N_B)` exactly.

3.  Max-statistic over the null controls FWER in the strong sense.

Method-specific step:

- `method = "threshold"` (default): threshold at
  `|t| > cluster_threshold`, label connected components with
  4-connectivity, cluster mass = sum of t-values within the cluster.

- `method = "tfce"`: enhance the observed t-map (and each permutation
  t-map) into a TFCE map; record `max(|TFCE|)` per permutation.
  Per-pixel p-value = fraction of null max-TFCE that equals or exceeds
  the observed `|TFCE|`. No free threshold parameter.

## Reading the result (threshold method)

- `$observed_t`, per-pixel Welch t vector.

- `$clusters`, data.frame with `cluster_id`, `direction`, `mass`,
  `size`, `p_value`, `significant`.

- `$null_distribution$pos`, `$neg`, per-permutation max masses.

- `$pos_labels`, `$neg_labels`, integer matrices for plotting.

- `$cluster_threshold`, `$alpha`, `$n_permutations`,
  `$n_participants_a`, `$n_participants_b`, `$method = "threshold"`.

## Reading the result (TFCE method)

- `$observed_t`, per-pixel Welch t vector (before enhancement).

- `$tfce_map`, per-pixel signed TFCE values.

- `$tfce_pmap`, per-pixel p-values against the max-TFCE null.

- `$tfce_significant_mask`, logical vector flagging pixels with
  `p < alpha`.

- `$null_distribution$max_abs_tfce`, per-permutation max\|TFCE\|.

- `$tfce_H`, `$tfce_E`, `$tfce_n_steps`, `$alpha`, `$n_permutations`,
  `$n_participants_a`, `$n_participants_b`, `$method = "tfce"`.

## Reading the plot

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the returned
object renders the per-pixel observed Welch t (or signed TFCE values for
`method = "tfce"`) as a raster image with FWE-significant cluster
boundaries traced on top.

- **Colour** encodes the sign and magnitude of the per-pixel statistic.
  Blue = condition A larger than B at that pixel (positive Welch t); red
  = condition B larger than A (negative Welch t); white = near zero. The
  colourbar reads in `Welch t` units (or signed `TFCE value` for TFCE).

- **Black contours** trace the boundary of clusters that are significant
  at the chosen `alpha` under the max-mass null (threshold method) or
  the boundary of pixels with FWE- corrected `p < alpha` (TFCE method).
  Pixels inside the contour survive the FWE correction; pixels outside
  do not.

- The colour scale is symmetric around zero by default so the neutral
  colour aligns with `t = 0`. The displayed t-map is the *observed*
  statistic before any thresholding; the contours encode the inferential
  decision.

- Colour convention matches
  [`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md)
  and
  [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)
  so the same data reads consistently across the package.

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) accepts
  `main =` (title) and `colour_bar =` (logical, default `TRUE`). Set
  `colour_bar = FALSE` when packing many panels into a small grid — used
  internally by
  [`plot.rcisignal_rel_pairwise_report()`](https://olivethree.github.io/rcisignal/reference/plot.rcisignal_rel_pairwise_report.md)
  for that reason.

## Reliability metrics expect raw masks

Welch t and cluster mass / TFCE are variance-based and sensitive to any
scaling. Inputs with `attr(., "source") == "rendered"` (set
automatically by Mode 1 readers like
[`extract_signal()`](https://olivethree.github.io/rcisignal/reference/extract_signal.md))
error unless `acknowledge_scaling = TRUE` cascades through the internal
[`pixel_t_test()`](https://olivethree.github.io/rcisignal/reference/pixel_t_test.md)
call.

## References

Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of
EEG- and MEG-data. *Journal of Neuroscience Methods*, 164(1), 177-190.
[doi:10.1016/j.jneumeth.2007.03.024](https://doi.org/10.1016/j.jneumeth.2007.03.024)

Nichols, T. E., & Holmes, A. P. (2002). Nonparametric permutation tests
for functional neuroimaging: a primer with examples. *Human Brain
Mapping*, 15(1), 1-25.
[doi:10.1002/hbm.1058](https://doi.org/10.1002/hbm.1058)

Smith, S. M., & Nichols, T. E. (2009). Threshold-free cluster
enhancement: addressing problems of smoothing, threshold dependence and
localisation in cluster inference. *NeuroImage*, 44(1), 83-98.
[doi:10.1016/j.neuroimage.2008.03.061](https://doi.org/10.1016/j.neuroimage.2008.03.061)

## See also

[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md),
[`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal call-signature demo with two synthetic inputs.
n_side <- 32L
n_pix  <- n_side * n_side
n_prod <- 20L
set.seed(1)
signal_matrix_a <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
signal_matrix_b <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
rel_cluster_test(signal_matrix_a, signal_matrix_b,
                 img_dims = c(n_side, n_side),
                 n_permutations = 200L, seed = 1)
} # }

if (FALSE) { # \dontrun{
# Same function, richer input: signal planted in different face regions
# (eyes vs mouth). The cluster test should find significant clusters
# at both regions, with opposite signs.
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
ct <- rel_cluster_test(sig_eyes, sig_mouth,
                       n_permutations = 500L, seed = 1)
# Plot the t-map with significant cluster boundaries: blue = eyes
# condition larger; red = mouth condition larger; black contours =
# FWE-significant clusters at alpha.
plot(ct, main = "Eyes vs Mouth (cluster-FWE)")
} # }
```
