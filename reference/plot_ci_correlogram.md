# Pairwise Pearson-r correlogram across multiple group-mean CIs

Renders a publication-ready Pearson correlation matrix across two or
more classification images (CIs), with optional face region masking and
optional direct save to PNG (600 dpi) or PDF.

Designed for the common publication-figure task: given a small set of
group-mean CIs from different conditions, show how similar they are
pixel-wise. The function accepts CIs in three shapes (column-major
vectors, single-column matrices, or per-producer `signal_matrix`
returned by `ci_from_responses_*()`); per-producer matrices are reduced
to group means via [`rowMeans()`](https://rdrr.io/r/base/colSums.html)
automatically.

Diverging-palette convention matches the rest of the package: positive
`r` = blue, negative `r` = red, neutral colour at zero. The colour scale
is fixed at `c(-1, 1)` so panels are directly comparable across runs and
across paper figures.

## Usage

``` r
plot_ci_correlogram(
  cis,
  img_dims = NULL,
  mask = c("none", "face", "upper_face", "lower_face"),
  triangle = c("full", "upper", "lower"),
  palette = c("diverging", "diverging_puor", "diverging_brbg"),
  show_values = TRUE,
  value_digits = 2L,
  main = NULL,
  file = NULL,
  width = NULL,
  height = NULL,
  ...
)
```

## Arguments

- cis:

  A named list of CIs to correlate. Each element may be a numeric vector
  of length `prod(img_dims)` (a single group-mean CI in column-major
  order), a numeric matrix `n_pixels x 1`, or a numeric matrix
  `n_pixels x n_producers` (a per-producer `signal_matrix` from
  `ci_from_responses_*()`; the function takes
  [`rowMeans()`](https://rdrr.io/r/base/colSums.html) to obtain the
  group mean). Names become the row / column / diagonal labels in the
  figure. At least two elements are required.

- img_dims:

  Integer `c(nrow, ncol)`. If `NULL`, inferred from
  `attr(cis[[1]], "img_dims")` (set by `ci_from_responses_*()`) or from
  `sqrt(n_pixels)` if that is a whole number.

- mask:

  One of `"none"` (default; correlate over all pixels), `"face"`
  (full-face oval from
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)),
  `"upper_face"` (top half of the oval), or `"lower_face"` (bottom
  half). Restricting to the face oval removes off-face noise from the
  correlation; sub-regions answer "do these CIs covary in the upper /
  lower face only".

- triangle:

  One of `"full"` (default; render the whole matrix), `"upper"` (mask
  the lower triangle and diagonal; put CI labels on the diagonal cells),
  or `"lower"` (mirror).

- palette:

  One of `"diverging"` (default; RdBu, blue = positive),
  `"diverging_puor"` (PuOr, purple = positive), or `"diverging_brbg"`
  (BrBG, green = positive). All three are colorblind-friendly diverging
  palettes appropriate for correlation matrices.

- show_values:

  Logical. If `TRUE` (default), render the correlation values inside
  each visible cell.

- value_digits:

  Integer. Decimal places for the cell-value labels. Default `2L` (e.g.,
  `"+0.67"`).

- main:

  Optional plot title.

- file:

  Optional output path. If `NULL` (default), plots to the current open
  device. If a path ending in `.png` or `.pdf` (case-insensitive), opens
  the corresponding device and writes the figure to disk: PNG at 600 dpi
  (raster), PDF as vector. Anything else aborts.

- width, height:

  Optional output dimensions in inches. Default sizes the canvas to be
  square and scale with the number of CIs:
  `max(5.5, 0.55 * n_cis + 4.5)` inches per side. Override either to
  control the final figure size.

- ...:

  Currently unused; reserved for future arguments.

## Value

Invisibly, a list with `correlation_matrix` (the n x n named Pearson-r
matrix actually plotted), `n_pixels_used`, `mask`, `palette`,
`triangle`, and `file` (the path written to, or `NULL`).

## Reading the plot

Each cell shows the Pearson correlation `r` between two group-mean CIs
computed across the pixels included in `mask` (all pixels by default).
Saturation encodes `|r|`; hue encodes sign (blue = positive, red =
negative). Numbers in each cell show the exact `r` to `value_digits`
decimal places.

The diagonal is identically 1 by definition. With `triangle = "full"`
the diagonal is rendered as a saturated blue cell with `+1.00` inside.
With `triangle = "upper"` or `triangle = "lower"` the diagonal cells are
blank and carry the CI labels instead, freeing the figure margins.

## Pearson r between CIs is not a clean similarity score

Two base-subtracted CIs share image-domain spatial structure (face
shape, oval signal support) that pushes their correlation above zero
even when the underlying mental representations are unrelated. A high
absolute `r` does not by itself license a similarity claim. Use this
correlogram to visualise the **relative** ordering of pair-wise
correlations (which pairs covary more than others), and benchmark any
absolute claim against a permutation null. For a baseline-free magnitude
summary of how different two conditions are, see
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
(Euclidean distance with bootstrap CI). For pixel-level localisation of
*where* two conditions differ, see
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md).

## See also

[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
for a Euclidean-distance (baseline-free) magnitude summary between two
conditions;
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)
for pixel-level inference on where two conditions differ;
[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
for the masks used by `mask = `.

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal call-signature demo with synthetic vectors.
n_side <- 32L
n_pix  <- n_side * n_side
set.seed(1)
ci_list <- list(
  A = rnorm(n_pix),
  B = rnorm(n_pix),
  C = rnorm(n_pix)
)
plot_ci_correlogram(ci_list, img_dims = c(n_side, n_side))
} # }

if (FALSE) { # \dontrun{
# Realistic input: simulate three conditions with planted signals
# in different face regions, then compare their CIs.
sim_eyes  <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "eyes",
  signal_region = "eyes", signal_strength = "strong", seed = 1
)
sim_mouth <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "mouth",
  signal_region = "mouth", signal_strength = "strong", seed = 2
)
sim_nose  <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "nose",
  signal_region = "nose", signal_strength = "strong", seed = 3
)
cis_eyes  <- ci_from_responses_briefrc(sim_eyes$data,
                                       noise_matrix = sim_eyes$noise_matrix)
cis_mouth <- ci_from_responses_briefrc(sim_mouth$data,
                                       noise_matrix = sim_mouth$noise_matrix)
cis_nose  <- ci_from_responses_briefrc(sim_nose$data,
                                       noise_matrix = sim_nose$noise_matrix)

plot_ci_correlogram(
  list("Eyes" = cis_eyes$signal_matrix,
       "Mouth" = cis_mouth$signal_matrix,
       "Nose"  = cis_nose$signal_matrix),
  mask     = "face",
  triangle = "upper",
  file     = "ci_correlogram.pdf"
)
} # }
```
