# Pairwise Euclidean distance matrix across multiple group-mean CIs

Renders a publication-ready Euclidean distance matrix across two or more
classification images (CIs). A stage-2 consumer: accepts a
[`group_ci()`](https://olivethree.github.io/rcisignal/reference/group_ci.md)
result directly, a named numeric matrix of `pixels x n_groups`, or the
same named-list-of-CIs format as
[`plot_ci_correlogram()`](https://olivethree.github.io/rcisignal/reference/plot_ci_correlogram.md)
(per-producer matrices reduced to group means internally).

Use this when the question is "how *far apart* are these CIs in pixel
space?" rather than "how do they covary?". Euclidean distance does not
share the positive-baseline issue Pearson `r` has on base-subtracted CIs
(see
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
for the two-condition bootstrap version of the same magnitude metric).

## Usage

``` r
plot_ci_distance_matrix(
  cis,
  img_dims = NULL,
  mask = c("none", "face", "upper_face", "lower_face"),
  method = c("raw", "normalised"),
  triangle = c("full", "upper", "lower"),
  palette = c("viridis", "inferno", "plasma", "rocket"),
  show_values = TRUE,
  value_digits = NULL,
  main = NULL,
  file = NULL,
  width = NULL,
  height = NULL,
  ...
)
```

## Arguments

- cis:

  CIs to compare. The recommended form is a numeric matrix
  `n_pixels x n_cis` with named columns; each column is one CI (a group
  mean, a single producer's CI, or any mix). Build it outside the call
  with `cbind(name = rowMeans(cis$signal_matrix), ...)`, or use the
  output of
  [`group_ci()`](https://olivethree.github.io/rcisignal/reference/group_ci.md)
  (a named matrix) directly. A named list of CIs is also accepted
  (vectors of length `prod(img_dims)`, single-column matrices, or
  per-producer `signal_matrix` objects from `ci_from_responses_*()`,
  which are reduced to group means internally). Names become the row /
  column / diagonal labels in the figure. At least two CIs are required.

- img_dims:

  Integer `c(nrow, ncol)`. If `NULL`, inferred from
  `attr(cis[[1]], "img_dims")` or from `sqrt(n_pixels)`.

- mask:

  One of `"none"` (default), `"face"`, `"upper_face"`, or
  `"lower_face"`. Restricts the pixels included in the distance
  computation via
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md).

- method:

  One of `"raw"` (default; absolute Euclidean) or `"normalised"`
  (`raw / sqrt(n_pixels_used)`).

- triangle:

  One of `"full"` (default), `"upper"` (mask the lower triangle and
  diagonal; CI labels on the diagonal cells), or `"lower"`.

- palette:

  One of `"viridis"` (default), `"inferno"`, `"plasma"`, or `"rocket"`.
  All four are sequential, colorblind-friendly palettes appropriate for
  non-negative distances.

- show_values:

  Logical. If `TRUE` (default), render the distance values inside each
  visible cell.

- value_digits:

  Integer or `NULL`. Decimal places for the cell labels. When `NULL`
  (default): 2 for raw distance, 3 for normalised distance.

- main:

  Optional plot title.

- file:

  Optional output path. If `NULL` (default), plots to the current open
  device. If a path ending in `.png` or `.pdf` (case-insensitive), saves
  at 600 dpi (PNG) or as vector PDF.

- width, height:

  Optional output dimensions in inches. Defaults size the canvas to be
  square and scale with the number of CIs.

- ...:

  Currently unused; reserved for future arguments.

## Value

Invisibly, a list with `distance_matrix`, `distance_raw`, `method`,
`n_pixels_used`, `mask`, `palette`, `triangle`, and `file`. See "Reading
the result".

## Reading the plot

Each cell shows the Euclidean distance between two group-mean CIs
computed across the pixels included in `mask` (all pixels by default).
Larger values mean the two CIs sit farther apart in pixel space.
Saturation encodes magnitude on a single-hue sequential ramp: pale (or
near-white) at distance ~ 0, deep / dark at large distances. Numbers in
each cell show the exact distance to `value_digits` decimal places.

The diagonal is identically 0 (a CI's distance from itself). With
`triangle = "full"` it is rendered as the lightest cell. With
`triangle = "upper"` or `triangle = "lower"` the diagonal cells are
blank and carry the CI labels instead, freeing the figure margins for
paper figures.

## Raw vs normalised distance

`method = "raw"` (default) returns the absolute Euclidean distance, in
whatever units the CIs carry. Useful for comparing contrasts computed on
the same image at the same resolution.

`method = "normalised"` divides by `sqrt(n_pixels_used)`, producing a
per-pixel root-mean-square distance that is comparable across masks and
resolutions. Use this when the question is "is contrast A more separated
than contrast B" and the two contrasts were computed under different
masks (or you want a publication number whose order of magnitude does
not balloon with image size).

Both matrices are always returned (see "Reading the result" below); the
`method` argument only controls which is rendered in the figure and used
for cell labels.

## Reading the result

The function returns invisibly a list with:

- `$distance_matrix`: the n x n distance matrix actually rendered (raw
  or normalised per `method`). Rows and columns are named with the CI
  names.

- `$distance_raw`: the raw distance matrix (always; useful for
  downstream MDS / [`hclust()`](https://rdrr.io/r/stats/hclust.html)
  regardless of which the figure showed).

- `$method`, `$n_pixels_used`, `$mask`, `$palette`, `$triangle`,
  `$file`: bookkeeping.

## See also

[`plot_ci_correlogram()`](https://olivethree.github.io/rcisignal/reference/plot_ci_correlogram.md)
for the Pearson-r version of the same input format;
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
for a two-condition Euclidean distance with bootstrap CI;
[`plot_ci_mds()`](https://olivethree.github.io/rcisignal/reference/plot_ci_mds.md)
to project multiple CIs into a 2D MDS scatter using the same distance
matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal: synthetic CIs to see the function's call signature
# and inspect the output shape.
set.seed(1)
n_pix <- 32L * 32L
ci_list <- list(
  A = rnorm(n_pix),
  B = rnorm(n_pix) + 0.3,
  C = rnorm(n_pix) - 0.2,
  D = rnorm(n_pix) + 0.5
)
out <- plot_ci_distance_matrix(ci_list, img_dims = c(32L, 32L))

# The distance matrix (named numeric matrix):
out$distance_matrix

# Raw distance is always returned, even if `method = "normalised"`:
out$distance_raw
} # }

if (FALSE) { # \dontrun{
# Realistic: simulate three conditions with planted signals in
# different face regions, build CIs, then compare in distance space.
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

plot_ci_distance_matrix(
  list("Eyes" = cis_eyes$signal_matrix,
       "Mouth" = cis_mouth$signal_matrix,
       "Nose"  = cis_nose$signal_matrix),
  mask     = "face",
  method   = "normalised",   # comparable across masks / resolutions
  triangle = "upper",
  file     = "ci_distance_matrix.pdf"
)
} # }
```
