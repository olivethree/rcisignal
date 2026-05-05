# Compare Euclidean dissimilarity intervals across multiple contrasts on a single plot

Side-by-side comparison of
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
results: one point per contrast, with its 95% percentile bootstrap CI as
a horizontal range. Useful for paper figures showing whether two
contrasts (e.g., Trust vs Friendly, Dominant vs Competent) have
overlapping or non-overlapping CIs without forcing the reader to read
four numbers from a table.

## Usage

``` r
plot_dissimilarity_grid(
  ...,
  metric = c("euclidean", "euclidean_normalised"),
  main = "Between-condition Euclidean distance"
)
```

## Arguments

- ...:

  Named `rcisignal_rel_dissim` objects (use `name = obj` syntax). Names
  become the y-axis labels.

- metric:

  `"euclidean"` (default) or `"euclidean_normalised"`.

- main:

  Plot title.

## Value

Invisibly a data frame of (label, observed, ci_low, ci_high) for further
use.

## See also

[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
n_pix  <- 32L * 32L
n_prod <- 20L
set.seed(1)
sig_a <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
sig_b <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
sig_c <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)

d_AB <- rel_dissimilarity(sig_a, sig_b, n_boot = 200L, seed = 1)
d_AC <- rel_dissimilarity(sig_a, sig_c, n_boot = 200L, seed = 1)
plot_dissimilarity_grid("A vs B" = d_AB, "A vs C" = d_AC)
} # }
```
