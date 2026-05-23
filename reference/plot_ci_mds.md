# Multidimensional-scaling (MDS) projection of multiple CIs

Projects two or more classification images (CIs) into a low- dimensional
Euclidean scatter using classical MDS
([`stats::cmdscale`](https://rdrr.io/r/stats/cmdscale.html)). Distances
between points in the scatter reproduce the original Euclidean distances
between CIs as faithfully as the chosen number of dimensions allows.

Use this to see at a glance which condition CIs cluster together and
which sit far apart, especially when you have more than a handful of
conditions to compare side by side.

## Usage

``` r
plot_ci_mds(
  cis,
  img_dims = NULL,
  mask = c("none", "face", "upper_face", "lower_face"),
  distance = c("euclidean_raw", "euclidean_normalised"),
  k = "auto",
  stress_threshold = 0.05,
  k_max = 4L,
  groups = NULL,
  shapes = NULL,
  point_size = 2.5,
  point_alpha = 0.85,
  label = TRUE,
  label_cex = 0.85,
  show_axes_at_zero = TRUE,
  show_gof = TRUE,
  main = NULL,
  file = NULL,
  width = NULL,
  height = NULL,
  ...
)
```

## Arguments

- cis:

  CIs to project, in any of three forms:

  - a
    [`group_ci()`](https://olivethree.github.io/rcisignal/reference/group_ci.md)
    result (an
    [rcisignal_group_ci](https://olivethree.github.io/rcisignal/reference/group_ci.md)
    matrix);

  - a numeric matrix `n_pixels x n_groups` with named columns;

  - a named list of CIs (each element a vector of length
    `prod(img_dims)`, a single-column matrix, or a per-producer
    `signal_matrix` from `ci_from_responses_*()` which is reduced to a
    group mean internally). At least three CIs are required for a
    meaningful 2D projection.

- img_dims:

  Integer `c(nrow, ncol)`. If `NULL`, inferred from
  `attr(cis[[1]], "img_dims")` or from `sqrt(n_pixels)`.

- mask:

  One of `"none"` (default), `"face"`, `"upper_face"`, or
  `"lower_face"`. Restricts the pixels included in the distance
  computation via
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md).

- distance:

  One of `"euclidean_raw"` (default; absolute Euclidean distance) or
  `"euclidean_normalised"` (`raw / sqrt(n_pixels_used)`). The MDS
  projection is fit to whichever distance is selected.

- k:

  Either `"auto"` (default; pick the smallest `k` whose Kruskal stress
  reaches `stress_threshold`) or an integer (forces that many
  dimensions; bypasses auto- selection). Use `k = 2L` for a single
  paper-figure panel.

- stress_threshold:

  Numeric. The Kruskal stress level the auto-selector treats as "good
  enough". Default `0.05` (Kruskal's "good" band). Lower it to `0.025`
  for "excellent" fidelity; raise to `0.10` for "fair" if the data are
  too noisy to reach 0.05 at small `k`.

- k_max:

  Integer. Maximum dimensionality the auto-selector will try. Default
  `4L`, capping the multi-panel grid at `choose(4, 2) = 6` pairs.
  Silently capped to `length(cis) - 1L` (the maximum embeddable
  dimension).

- groups:

  Optional named character vector mapping CI names to a categorical
  group label. Used to colour the points. Names must match `names(cis)`
  exactly.

- shapes:

  Optional named character vector with the same shape as `groups`,
  mapping CI names to a second categorical label used for point shapes
  (`pch`).

- point_size, point_alpha:

  Numeric. Point cex and alpha.

- label:

  Logical. Draw CI names above each point. Default `TRUE`.

- label_cex:

  Numeric. Cex for point labels.

- show_axes_at_zero:

  Logical. Draw faint zero-reference lines on each panel. Default
  `TRUE`.

- show_gof:

  Logical. Render the GOF summary above the figure (selected `k`,
  stress + band, per-axis variance). Default `TRUE`.

- main:

  Optional plot title.

- file:

  Optional output path. If `NULL` (default), plots to the current open
  device. If a path ending in `.png` or `.pdf` (case-insensitive), saves
  at 600 dpi (PNG) or as vector PDF.

- width, height:

  Optional output dimensions in inches. Defaults scale with the number
  of CIs and the number of panels.

- ...:

  Currently unused; reserved for future arguments.

## Value

Invisibly, an object of class `rcisignal_mds`. See "Reading the result".

## Auto-selecting the number of dimensions

By default (`k = "auto"`), the function fits classical MDS at every
dimensionality from 2 up to `k_max` and picks the **smallest** `k` whose
Kruskal stress-1 against the original Euclidean distances reaches the
"good" threshold (default `stress_threshold = 0.05`).

Kruskal's (1964) interpretive bands for stress-1:

- `0.025` excellent (the 2D / kD map is essentially exact)

- `0.05` good (small distortions; safe to interpret)

- `0.10` fair (interpret carefully; check the trace)

- `0.20` poor (the projection is hiding more than it shows)

- `> 0.20` very poor

When the selected `k > 2`, the figure becomes a grid of all
`choose(k, 2)` pairwise dimension panels (Dim 1 vs 2, Dim 1 vs 3, Dim 2
vs 3, ...) so no information is hidden by a premature flattening to 2D.

The companion goodness-of-fit metric (Mardia's GOF1, cumulative
eigenvalue variance) is reported alongside stress in
`$variance_pct_by_k`. Values near 100% mean the kD map captures
essentially all of the positive eigenmass; below ~60% the kD map is
hiding more than it shows. Stress is usually the more interpretable of
the two.

## Theory-driven dimensionality

Real research often has a theoretical reason to fix `k` to a specific
number (typically 2, for a paper-style scatter that matches a two-axis
hypothesis like warmth/dominance, or a single panel for visual
comparison with prior work). Pass an integer `k` to bypass
auto-selection:

    plot_ci_mds(ci_list, mask = "face", k = 2L)

The function still computes stress at every `k` in `[2, k_max]` and
exposes the trace via `$stress_by_k`. When the forced dimensionality has
high stress, report `out$stress_by_k` alongside the 2D figure in the
paper so readers see what the theory-driven projection is hiding;
interpret point positions in *relative* terms (which conditions cluster
together) rather than as absolute distances.

Force a higher specific `k` (e.g., `k = 3L`) when a richer projection is
theoretically motivated; the function will render the corresponding
`choose(k, 2)`-panel grid.

## Reading the plot

Each panel is a scatter of CIs in two MDS dimensions. Points close
together represent CIs that are similar in pixel space; distant points
represent dissimilar CIs. With `groups`, points colour by category; with
`shapes`, points carry distinct marker shapes by another category. The
aspect ratio is fixed (`asp = 1`) so visual distances faithfully reflect
the MDS distances.

The line at the top of the figure reports the selected `k`, the Kruskal
stress at that `k` and its interpretive band, and the per-axis variance
percentages.

## Reading the result

The returned object (class `rcisignal_mds`) contains:

- `$mds_points`: **the coordinates of each CI in the Euclidean MDS
  space**. An `n_cis x k_selected` named numeric matrix (rows = CI
  names, columns = orthogonal MDS axes). Pull these directly to compute
  custom downstream analyses (e.g. distances between specific CIs in the
  reduced space, your own scatter, clustering).

- `$distance_matrix`: the pairwise distances the projection was fit to
  (raw or normalised per `distance`).

- `$stress_by_k`, `$variance_pct_by_k`: per-`k` Kruskal stress and
  cumulative Mardia GOF1. Inspect these to audit the auto-selected `k`.

- `$k_selected`, `$stress_1`, `$stress_band`, `$variance_pct`: summaries
  at the selected `k`.

- `$panel_pairs`: the `2 x choose(k, 2)` matrix of dimension pairs the
  figure rendered.

Call [`print()`](https://rdrr.io/r/base/print.html) or
[`summary()`](https://rdrr.io/r/base/summary.html) on the returned
object for a one-screen human-readable view.

## References

Kruskal, J. B. (1964). Multidimensional scaling by optimizing goodness
of fit to a nonmetric hypothesis. *Psychometrika*, 29(1), 1-27.

## See also

[`plot_ci_distance_matrix()`](https://olivethree.github.io/rcisignal/reference/plot_ci_distance_matrix.md)
for the all-vs-all distance matrix that this function projects;
[`plot_ci_correlogram()`](https://olivethree.github.io/rcisignal/reference/plot_ci_correlogram.md)
for the Pearson-r version of the same multi-CI summary;
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
for a two-condition Euclidean distance with bootstrap CI.

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
  D = rnorm(n_pix) + 0.5,
  E = rnorm(n_pix)
)

# Auto-selects the smallest k reaching "good" stress (<= 0.05).
out <- plot_ci_mds(ci_list, img_dims = c(32L, 32L))

# Human-readable view of the dimensionality selection trace:
print(out)

# Extract the coordinates of each CI in the Euclidean MDS space.
out$mds_points

# Distance between CI "A" and CI "C" in the reduced space:
sqrt(sum((out$mds_points["A", ] - out$mds_points["C", ])^2))

# Per-k stress and variance traces (was the auto-choice sensible?):
out$stress_by_k
out$variance_pct_by_k
} # }

if (FALSE) { # \dontrun{
# Realistic: simulate four conditions with planted signals,
# build CIs, then compare in MDS space.
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
sim_flat  <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "control",
  signal_region = NULL, seed = 4
)

cis_eyes  <- ci_from_responses_briefrc(sim_eyes$data,
                                       noise_matrix = sim_eyes$noise_matrix)
cis_mouth <- ci_from_responses_briefrc(sim_mouth$data,
                                       noise_matrix = sim_mouth$noise_matrix)
cis_nose  <- ci_from_responses_briefrc(sim_nose$data,
                                       noise_matrix = sim_nose$noise_matrix)
cis_flat  <- ci_from_responses_briefrc(sim_flat$data,
                                       noise_matrix = sim_flat$noise_matrix)

ci_list <- list(
  "Eyes"    = cis_eyes$signal_matrix,
  "Mouth"   = cis_mouth$signal_matrix,
  "Nose"    = cis_nose$signal_matrix,
  "Control" = cis_flat$signal_matrix
)

out <- plot_ci_mds(ci_list, mask = "face")

# Force a single 2D paper figure once you have audited fidelity:
plot_ci_mds(ci_list, mask = "face", k = 2L,
            file = "fig_mds.pdf")

# Add a grouping variable for point colour:
plot_ci_mds(
  ci_list, mask = "face",
  groups = c(Eyes = "feature", Mouth = "feature",
             Nose = "feature", Control = "control")
)
} # }
```
