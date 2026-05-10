# Run all pairwise between-condition comparisons across K conditions

Generalises
[`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md)
from a 2-condition comparison to K conditions: runs
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)
and
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
on every K-choose-2 pair and applies a family-wise error correction
across pairs.

## Usage

``` r
run_discriminability_pairwise(
  signal_matrices,
  fwer = c("holm", "bonferroni", "none"),
  img_dims = NULL,
  paired = FALSE,
  method = c("threshold", "tfce"),
  n_permutations = 2000L,
  n_boot = 2000L,
  cluster_threshold = 2,
  alpha = 0.05,
  ci_level = 0.95,
  dissim_null = c("none", "permutation"),
  mask = NULL,
  seed = NULL,
  progress = TRUE,
  acknowledge_scaling = FALSE
)
```

## Arguments

- signal_matrices:

  Named list of pixels x participants signal matrices, one per
  condition. Names become condition labels in the output.

- fwer:

  One of `"holm"` (default), `"bonferroni"`, or `"none"`.

- img_dims:

  Integer `c(nrow, ncol)`. If `NULL`, inferred from
  `attr(signal_matrices[[1]], "img_dims")`.

- paired:

  Logical. When `TRUE`, all pairs use the paired variant. All matrices
  must have identical column names.

- method:

  Cluster-test method. Default `"threshold"`.

- n_permutations, n_boot, cluster_threshold, alpha, ci_level:

  Forwarded to per-pair calls.

- dissim_null:

  Forwarded to
  [`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
  as `null`. Default `"none"` to keep wall time bounded.

- mask:

  Optional logical vector of length `nrow(signal_matrices[[1]])`
  (column-major) forwarded to both per-pair calls. Build with
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  (parametric oval and sub-regions) or
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
  (PNG/JPEG mask).

- seed:

  Optional integer.

- progress:

  Show `cli` progress bars per per-pair call.

- acknowledge_scaling:

  Logical. Forwarded.

## Value

Object of class `rcisignal_rel_pairwise_report` with:

- `$pairs`, data.frame: `pair_id`, `cond_a`, `cond_b`, `n_clusters`,
  `p_min`, `p_adj_pair`, `significant`, `euclidean`,
  `euclidean_normalised`.

- `$results`, named list of per-pair `rel_cluster_test` and
  `rel_dissimilarity` results.

- `$conditions`, the input names.

- `$fwer`, the correction method used.

- `$alpha`, the across-pairs alpha.

## FWER scope

Cluster-level p-values within each pair are already max-statistic
FWER-controlled by
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md).
The Holm/Bonferroni adjustment here controls family-wise error across
the K-choose-2 pair comparisons (a second layer above the cluster test's
internal control), not over individual pixels or clusters.

For each pair, the statistic carried into the across-pairs adjustment is
the minimum cluster-level p-value within that pair. Within-pair cluster
p-values are not re-adjusted; they remain the max-statistic
FWER-controlled values from
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md).
A pair with no clusters contributes `p_min = 1.0` so the Holm ordering
is well-defined.

## Reading the plot

[`print()`](https://rdrr.io/r/base/print.html) shows the across-pairs
table (one row per pair, adjusted p-values, Euclidean magnitudes). To
visualise a single pair, pull its child results out of `$results` and
call [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on either:
e.g. `plot(res$results[["A_vs_B"]]$cluster_test)` for the cluster t-map
(blue = A larger; red = B larger; black contours bound FWE-significant
clusters), or `plot(res$results[["A_vs_B"]]$dissimilarity)` for the
bootstrap distribution. To compare Euclidean distances across pairs on a
shared axis, pass each pair's `dissimilarity` child to
[`plot_dissimilarity_grid()`](https://olivethree.github.io/rcisignal/reference/plot_dissimilarity_grid.md).

## See also

[`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md),
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md),
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Three conditions, each with signal planted in a different face region
# (eyes / mouth / nose). All three pairwise contrasts should detect a
# difference; the FWER correction is applied across the three pairs.
sim_eyes  <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "x",
  signal_region = "eyes",  signal_strength = "strong", seed = 1
)
sim_mouth <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "x",
  signal_region = "mouth", signal_strength = "strong", seed = 2
)
sim_nose  <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "x",
  signal_region = "nose",  signal_strength = "strong", seed = 3
)
sigs <- list(
  eyes  = ci_from_responses_briefrc(
    sim_eyes$data,  noise_matrix = sim_eyes$noise_matrix)$signal_matrix,
  mouth = ci_from_responses_briefrc(
    sim_mouth$data, noise_matrix = sim_mouth$noise_matrix)$signal_matrix,
  nose  = ci_from_responses_briefrc(
    sim_nose$data,  noise_matrix = sim_nose$noise_matrix)$signal_matrix
)
res <- run_discriminability_pairwise(sigs,
                                     n_permutations = 500L,
                                     n_boot = 500L, seed = 1)
print(res)  # adjusted p-values across the three pairs
# Cluster t-map for one pair (blue = first condition larger;
# red = second condition larger; black contours = significant clusters):
plot(res$results[["eyes_vs_mouth"]]$cluster_test,
     main = "Eyes vs Mouth (cluster-FWE)")
# Compare Euclidean distances across all three pairs on a shared axis:
plot_dissimilarity_grid(
  "eyes vs mouth" = res$results[["eyes_vs_mouth"]]$dissimilarity,
  "eyes vs nose"  = res$results[["eyes_vs_nose"]]$dissimilarity,
  "mouth vs nose" = res$results[["mouth_vs_nose"]]$dissimilarity
)
} # }
```
