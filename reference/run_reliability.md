# Run every within-condition reliability metric

Convenience orchestrator that runs the reliability metrics proper
(split-half with Spearman-Brown correction and ICC) on a single
condition's signal matrix and wraps the two results in an
`rcisignal_rel_report` for joint printing and plotting.

Use this when you want the full within-condition reliability report in
one call. Pass each metric individually if you need to tune arguments
per metric.

## Usage

``` r
run_reliability(
  signal_matrix,
  n_permutations = 2000L,
  null = c("none", "permutation", "random_responders"),
  noise_matrix = NULL,
  icc_variants = c("3_1", "3_k"),
  mask = NULL,
  seed = NULL,
  progress = TRUE,
  acknowledge_scaling = FALSE
)
```

## Arguments

- signal_matrix:

  Pixels x participants, base-subtracted.

- n_permutations:

  Passed to
  [`rel_split_half()`](https://olivethree.github.io/rcisignal/reference/rel_split_half.md).
  Default 2000.

- null:

  Forwarded to
  [`rel_split_half()`](https://olivethree.github.io/rcisignal/reference/rel_split_half.md).
  Default `"none"`.

- noise_matrix:

  Required when `null = "random_responders"`; forwarded to
  [`rel_split_half()`](https://olivethree.github.io/rcisignal/reference/rel_split_half.md).

- icc_variants:

  Passed to
  [`rel_icc()`](https://olivethree.github.io/rcisignal/reference/rel_icc.md).

- mask:

  Optional logical vector of length `nrow(signal_matrix)`. Threaded
  through to both metrics.

- seed:

  Optional integer; used for the split-half permutations.

- progress:

  Show `cli` progress bars.

- acknowledge_scaling:

  Logical. Forwarded to
  [`rel_icc()`](https://olivethree.github.io/rcisignal/reference/rel_icc.md).

## Value

Object of class `rcisignal_rel_report` with `$results` = named list of
two result objects (`split_half`, `icc`) and `$method = "reliability"`.

## What is included (and what is not)

Returns the two metrics that quantify the reliability of the group-level
classification image proper: split-half (a permutation-based estimate of
group-CI stability with Spearman-Brown projection to the full sample)
and ICC(3,\*) (the psychometric variance decomposition). These are
non-redundant.

Leave-one-out influence screening lives in
[`rel_loo()`](https://olivethree.github.io/rcisignal/reference/rel_loo.md)
and is **not** bundled here. Its output is an influence diagnostic, not
a reliability statistic, and mixing it into the reliability report
invites mis-reading `r_loo` values (which are near 1 by construction) as
reliability.

## Reading the result

`$results$split_half`, `$results$icc`, one result object each, with the
same fields as the standalone functions. `$method = "reliability"`.

## Reliability metrics expect raw masks

Both downstream metrics expect the raw mask. Inputs with
`attr(., "source") == "rendered"` (set automatically by Mode 1 readers
like
[`extract_signal()`](https://olivethree.github.io/rcisignal/reference/extract_signal.md))
error in
[`rel_icc()`](https://olivethree.github.io/rcisignal/reference/rel_icc.md)
unless `acknowledge_scaling = TRUE`.

## See also

[`rel_split_half()`](https://olivethree.github.io/rcisignal/reference/rel_split_half.md),
[`rel_icc()`](https://olivethree.github.io/rcisignal/reference/rel_icc.md),
[`rel_loo()`](https://olivethree.github.io/rcisignal/reference/rel_loo.md)
for the influence diagnostic;
[`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# End-to-end demo: simulate -> compute CIs -> assess reliability.
sim <- simulate_briefrc_data(n_per_condition = 10, n_trials = 60, seed = 1)
cis <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
rel <- run_reliability(cis$signal_matrix,
                       n_permutations = 200L, seed = 1)
print(rel)
} # }
```
