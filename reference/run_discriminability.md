# Run every between-condition discriminability metric

Convenience orchestrator that runs
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)
and
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
on two condition signal matrices and wraps both results in an
`rcisignal_rel_report` for joint printing / plotting.

Use this when you want both the spatial-pattern test and the overall
magnitude test in one call.

## Usage

``` r
run_discriminability(
  signal_matrix_a,
  signal_matrix_b,
  img_dims = NULL,
  n_permutations = 2000L,
  n_boot = 2000L,
  cluster_threshold = 2,
  alpha = 0.05,
  ci_level = 0.95,
  mask = NULL,
  seed = NULL,
  progress = TRUE,
  acknowledge_scaling = FALSE
)
```

## Arguments

- signal_matrix_a, signal_matrix_b:

  Pixels x participants, base-subtracted. Row counts must match.

- img_dims:

  Integer `c(nrow, ncol)`. If `NULL`, inferred from
  `attr(signal_matrix_a, "img_dims")`.

- n_permutations:

  Passed to
  [`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md).
  Default 2000.

- n_boot:

  Passed to
  [`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md).
  Default 2000.

- cluster_threshold:

  Passed to
  [`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md).
  Default 2.0.

- alpha:

  Passed to
  [`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md).
  Default 0.05.

- ci_level:

  Passed to
  [`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md).
  Default 0.95.

- mask:

  Optional logical vector of length `nrow(signal_matrix_a)`
  (column-major) threaded through to both downstream calls. Build with
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  (parametric oval and sub-regions) or
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
  (PNG/JPEG mask).

- seed:

  Optional integer.

- progress:

  Show `cli` progress bars.

- acknowledge_scaling:

  Logical. Forwarded.

## Value

Object of class `rcisignal_rel_report` with `$results` = named list of
two result objects (`cluster_test`, `dissimilarity`) and
`$method = "discriminability"`.

## Reading the result

`$results$cluster_test` and `$results$dissimilarity`, one result object
each, fields as in the standalone functions.
`$method = "discriminability"`.

## Reading the plot

Calling `plot(result)` lays out one panel per child result (cluster
t-map with FWE-significant contours; bootstrap dissimilarity
histograms). To plot one panel at a time, call
`plot(result$results$cluster_test)` or
`plot(result$results$dissimilarity)` directly. Colour convention on the
cluster panel matches the rest of the package: blue = condition A
larger; red = condition B larger; black contours bound FWE-significant
clusters.

## Reliability metrics expect raw masks

Both downstream metrics are scale-sensitive: the cluster test uses
variance-based Welch t, and Euclidean distance in
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
is sensitive to any scaling. Inputs with
`attr(., "source") == "rendered"` error unless
`acknowledge_scaling = TRUE`.

## See also

[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md),
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md),
[`run_reliability()`](https://olivethree.github.io/rcisignal/reference/run_reliability.md),
[`run_discriminability_pairwise()`](https://olivethree.github.io/rcisignal/reference/run_discriminability_pairwise.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-condition pipeline with signal planted in different face regions
# (eyes vs mouth). The orchestrator runs the cluster test and the
# bootstrap dissimilarity in one call; significant pixels should
# localise around the eye and mouth regions.
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
res <- run_discriminability(sig_eyes, sig_mouth,
                             n_permutations = 500L, n_boot = 500L,
                             seed = 1)
print(res)
# Whole report (cluster panel + dissimilarity panel side by side):
plot(res)
# Or each panel on its own:
plot(res$results$cluster_test,
     main = "Eyes vs Mouth — cluster test")
plot(res$results$dissimilarity,
     main = "Eyes vs Mouth — dissimilarity")
} # }
```
