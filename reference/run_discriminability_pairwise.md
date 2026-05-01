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

  Optional logical vector forwarded to both per-pair calls.

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

## See also

[`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md),
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md),
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
