# Collapse a per-producer signal matrix into per-group means

Stage-2 aggregator of the rcisignal pipeline. Collapses a per-producer
`signal_matrix` (pixels x n_producers, the object returned by
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md)
/
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md))
into a per-group matrix (pixels x n_groups) for use with the
distance-matrix, MDS, and correlogram plot functions.

## Usage

``` r
group_ci(signal_matrix, by, drop = TRUE)
```

## Arguments

- signal_matrix:

  Numeric matrix of pixels x n_producers, as returned by
  [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md)
  or
  [`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md)
  in their `$signal_matrix` field. The column count must match
  `length(by)` (or `length(by[[1]])` when `by` is a list).

- by:

  Either an atomic vector of length `ncol(signal_matrix)` (one group
  label per producer) or a named list of such vectors for factorial
  grouping (cell names are the levels joined by `"_"`). Coerced to
  factor; `NA` levels are dropped with a warning naming the count.

- drop:

  Logical. When `TRUE` (default), empty cells are dropped from the
  output. When `FALSE`, empty cells are present as `NA` columns with
  `n = 0L`.

## Value

A numeric matrix of pixels x n_groups, classed
`c("rcisignal_group_ci", "matrix", "array")`. Carries:

- column names = group labels;

- `attr(., "n")` = named integer vector of per-group producer counts;

- `attr(., "img_dims")` = inherited from
  `attr(signal_matrix, "img_dims")` if present.

## Details

The package has two stages. Stage 1 (per-producer `signal_matrix`) is
the only object accepted by reliability, discriminability, and infoVal
functions. Stage 2 (group-averaged matrix) is for plotting, RDM-style
comparison, and MDS. `group_ci()` is the stage-1-to-stage-2 transformer.

For each group, the corresponding output column is
`rowMeans(signal_matrix[, producers_in_group, drop = FALSE])`. When `by`
is a list, the cell label is the levels joined by `"_"` in the order the
list elements are given.

`group_ci()` does not accept (and will never accept) `trial_counts`,
`noise_matrix`, or `mask`. Anything that needs producer-level
information lives in stage 1: do the analysis first, then aggregate.

## See also

Upstream (stage 1):
[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md),
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md).
Downstream (stage 2):
[`plot_ci_distance_matrix()`](https://olivethree.github.io/rcisignal/reference/plot_ci_distance_matrix.md),
[`plot_ci_mds()`](https://olivethree.github.io/rcisignal/reference/plot_ci_mds.md),
[`plot_ci_correlogram()`](https://olivethree.github.io/rcisignal/reference/plot_ci_correlogram.md).

## Examples

``` r
set.seed(1)
n_pix  <- 32L * 32L
n_prod <- 12L
sm <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod,
             dimnames = list(NULL, sprintf("p%02d", seq_len(n_prod))))
attr(sm, "img_dims") <- c(32L, 32L)
g <- rep(c("A", "B"), each = n_prod / 2L)
gcis <- group_ci(sm, by = g)
dim(gcis)                 # n_pix x 2
#> [1] 1024    2
attr(gcis, "n")           # named per-group producer counts
#> A B 
#> 6 6 

if (FALSE) { # \dontrun{
# Realistic end-to-end: simulate, build per-producer CIs, collapse.
sim <- simulate_briefrc_data(
  n_per_condition = 10, n_trials = 60,
  conditions = c("A", "B"), seed = 1
)
cis <- ci_from_responses_briefrc(sim$data,
                                 noise_matrix = sim$noise_matrix)
producer_to_cond <- sim$data$condition[match(
  colnames(cis$signal_matrix), sim$data$participant_id)]
gcis <- group_ci(cis$signal_matrix, by = producer_to_cond)
plot_ci_distance_matrix(gcis,
  img_dims = c(sim$meta$img_size, sim$meta$img_size))
} # }
```
