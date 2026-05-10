# Per-pixel agreement heatmap for a producer signal matrix

Visualises where producers in a single condition agree on the direction
of signal. For each pixel, computes a one-sample t-statistic against
zero across producers (`mean / (sd / sqrt(N))`), then displays the
resulting map with a diverging colour palette (positive = agreement on
positive signal, negative = agreement on negative signal, zero = no
agreement). Saturation of the colour is the magnitude of the agreement,
*not* the value of the group-mean CI.

Use this to answer "where do producers consistently *agree* the target
trait lives in the face?". Pair with the group-mean CI image (raw mask
or rendered) to see *direction* and *agreement* side by side.

## Usage

``` r
plot_agreement_map(
  signal_matrix,
  img_dims = NULL,
  mask = NULL,
  threshold = NULL,
  zlim = NULL,
  palette = c("diverging", "viridis"),
  main = "Per-pixel producer agreement (t-map)",
  ...
)
```

## Arguments

- signal_matrix:

  Pixels x participants raw mask (as returned by `ci_from_responses_*()`
  or
  [`read_cis()`](https://olivethree.github.io/rcisignal/reference/read_cis.md) +
  [`extract_signal()`](https://olivethree.github.io/rcisignal/reference/extract_signal.md)).

- img_dims:

  Integer `c(nrow, ncol)`. If `NULL`, inferred from
  `attr(signal_matrix, "img_dims")` or from `sqrt(n_pixels)` if the
  latter is a whole number.

- mask:

  Optional logical vector of length `nrow(signal_matrix)` (column-major)
  restricting display to a region (e.g.,
  `make_face_mask(img_dims, region = "eyes")`). Also accepts the output
  of
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
  for PNG/JPEG masks. Pixels outside the mask render as `NA`
  (transparent).

- threshold:

  Optional positive numeric. When supplied, pixels with
  `|t| < threshold` are rendered in the neutral (white) colour, making
  clusters of agreement stand out. Default `NULL` (full continuous map).

- zlim:

  Numeric `c(low, high)` for the colour scale. Default is symmetric
  around zero at `c(-max(|t|), max(|t|))` so the neutral colour aligns
  with t = 0.

- palette:

  Character. `"diverging"` (default; positive = blue, negative = red,
  neutral = white) or `"viridis"` (no neutral; for absolute-magnitude
  views).

- main:

  Title.

- ...:

  Passed to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html).

## Value

Invisibly, a list with `t_map` (numeric vector of t values per pixel),
`n` (producer count), `img_dims`, and `mask` (if supplied) — useful for
further analysis or replotting.

## Details

This is structurally a one-sample t-map (vs zero); pixels where
producers' contributions are large *and consistent in sign* get high
\|t\|, pixels where contributions are random get t near zero.
Cluster-permutation inference would normally accompany this for formal
pixel-level FWER control between conditions
([`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md));
the agreement map is the **descriptive** counterpart for a single
condition.

## Reading the plot

- **Colour** encodes the sign of the per-pixel one-sample t. Blue =
  producers consistently *add* to the base at that pixel (positive
  agreement); red = consistently *subtract* (negative agreement); white
  = no agreement (t near zero).

- **Saturation** encodes `|t|`: deeper colour means the agreement among
  producers is both large and consistent in sign. The colourbar on the
  right reads in `t` units (one-sample vs 0).

- **`zlim`** is symmetric around zero by default so the neutral colour
  aligns with `t = 0`. Pass `zlim = c(-z, z)` to fix the scale across
  panels for direct comparison.

- **`threshold`** clips colour to white below `|t| < threshold`, making
  strong-agreement clusters stand out. This is descriptive only; it does
  not provide FWER control. For inferential pixel significance, use
  [`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md)
  and overlay the contours via
  [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md).

- Colour convention matches
  [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)
  and the cluster-test plots so the same group CI reads consistently
  across the package.

## See also

[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md),
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)
for inferential between-condition tests.

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal call-signature demo with a synthetic input. The agreement
# map will look flat because the input is pure noise.
n_side <- 32L
n_pix  <- n_side * n_side
set.seed(1)
signal_matrix <- matrix(rnorm(n_pix * 20L), n_pix, 20L)
plot_agreement_map(signal_matrix, img_dims = c(n_side, n_side))
} # }

if (FALSE) { # \dontrun{
# Same function, richer input: simulate Brief-RC responses with a
# signal planted in the eye region, then look at the agreement map.
# Producers should consistently agree on the planted region.
sim <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "target",
  signal_region = "eyes", signal_strength = "strong", seed = 1
)
cis <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
plot_agreement_map(cis$signal_matrix)
} # }
```
