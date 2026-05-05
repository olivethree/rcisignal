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

  Optional logical vector of length `nrow(signal_matrix_a)`. Threaded
  through to both downstream calls.

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
# Two-condition pipeline: simulate -> per-condition CIs -> contrast.
sim <- simulate_briefrc_data(n_per_condition = 10, n_trials = 60,
                             conditions = c("target", "control"),
                             seed = 1)
a <- subset(sim$data, condition == "target")
b <- subset(sim$data, condition == "control")
sig_a <- ci_from_responses_briefrc(a, noise_matrix = sim$noise_matrix)$signal_matrix
sig_b <- ci_from_responses_briefrc(b, noise_matrix = sim$noise_matrix)$signal_matrix
run_discriminability(sig_a, sig_b,
                     n_permutations = 200L, n_boot = 200L, seed = 1)
} # }
```
