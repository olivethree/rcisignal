# Plot the cluster-test grid for a pairwise discriminability report

Renders one cluster t-map (or signed TFCE map) per pair, laid out in a
square-ish grid. Each panel shows the per-pixel statistic with black
contours bounding FWE-significant clusters (or FWE-corrected pixels
under TFCE). Colour convention matches the rest of the package: blue =
first condition larger, red = second condition larger.

## Usage

``` r
# S3 method for class 'rcisignal_rel_pairwise_report'
plot(x, ..., ncol = NULL, max_pairs = 12L, base_image = NULL, alpha_max = 0.7)
```

## Arguments

- x:

  A
  [`run_discriminability_pairwise()`](https://olivethree.github.io/rcisignal/reference/run_discriminability_pairwise.md)
  result.

- ...:

  Reserved for future use.

- ncol:

  Optional integer. Columns in the panel grid. When `NULL` (default),
  `ceiling(sqrt(n_pairs))` is used.

- max_pairs:

  Integer. Above this many pairs a warning is emitted (panels become
  illegible). Default `12L`. The grid is still drawn; pass
  `max_pairs = Inf` to silence the warning.

- base_image:

  Optional. Either a numeric matrix (`nrow x ncol`, grayscale, values in
  0-1) or a path to a PNG/JPEG file. When supplied, each per-pair t-map
  is composited on top of the grayscale base; otherwise panels render on
  a flat background.

- alpha_max:

  Numeric in `[0, 1]`. Maximum opacity of the t-map overlay at the
  colour-scale top when `base_image` is supplied. Default 0.7.

## Value

Invisibly the input `x`.

## Details

To compare overall magnitudes across pairs on a shared axis, pass the
per-pair `$dissimilarity` children to
[`plot_dissimilarity_grid()`](https://olivethree.github.io/rcisignal/reference/plot_dissimilarity_grid.md)
instead.

## See also

[`plot_dissimilarity_grid()`](https://olivethree.github.io/rcisignal/reference/plot_dissimilarity_grid.md)
for a shared-axis comparison of bootstrap dissimilarity distances across
pairs.
