# Per-pixel agreement heatmap for a producer signal matrix

Visualises where producers in a single condition agree on the direction
of signal. For each pixel, computes a one-sample t-statistic against
zero across producers (`mean / (sd / sqrt(N))`), then displays the
resulting map with a diverging color palette (positive = agreement on
positive signal, negative = agreement on negative signal, zero = no
agreement). Saturation of the color is the magnitude of the agreement,
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
  palette = c("diverging", "fire"),
  base_image = NULL,
  alpha_max = 0.7,
  main = "Per-pixel producer agreement (t-map)",
  show_n = TRUE,
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
  `|t| < threshold` are rendered in the neutral (white) color, making
  clusters of agreement stand out. Default `NULL` (full continuous map).

- zlim:

  Numeric `c(low, high)` for the color scale. For
  `palette = "diverging"` (default), defaults to
  `c(-max(|t|), max(|t|))` so the neutral color aligns with `t = 0`. For
  `palette = "fire"`, defaults to `c(0, max(|t|))` so pale yellow aligns
  with `|t| = 0`.

- palette:

  Character. `"diverging"` (default; positive = blue, negative = red,
  neutral = white) encodes sign in hue and magnitude in saturation.
  `"fire"` encodes `|t|` only on a single-hue ramp (pale yellow at zero
  -\> deep red at large `|t|`); use this when the question is "where do
  producers have a consistent opinion" and direction is not needed. The
  `"fire"` view discards sign; pair with `palette = "diverging"` or with
  [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)
  to recover direction at a region of interest.

- base_image:

  Optional. Either a numeric matrix (`nrow x ncol`, grayscale, values in
  0-1) or a path to a PNG/JPEG file. When supplied, the t-map is
  composited on top of the grayscale base; out-of-mask and subthreshold
  pixels render fully transparent. When `NULL` (default), the map is
  drawn on a flat panel via
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html) (the
  historical behavior).

- alpha_max:

  Numeric in `[0, 1]`. Maximum opacity of the heatmap at the color-scale
  top (`zlim_max`) when `base_image` is supplied. Ignored otherwise.
  Default 0.7.

- main:

  Title.

- show_n:

  Logical. When `TRUE` (default), draw the "N = ... producers, W x H
  pixels" subtitle line below the title. Set `FALSE` for multi-panel
  layouts where this information is already in the caption.

- ...:

  Passed to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html).

## Value

Invisibly, a list with `t_map` (numeric vector of t values per pixel;
always signed regardless of palette), `n` (producer count), `img_dims`,
`mask` (if supplied), `zlim` (the color scale used), and `palette` (the
palette name).

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

Two palettes are available; pick by what question you are asking the
data.

**`palette = "diverging"` (default).** Encodes sign and magnitude
together. Both deep red and deep blue indicate **strong** agreement
among producers; only the **direction** differs. "No agreement" is the
neutral color (white), not red.

- **Hue** encodes the sign of the per-pixel one-sample `t`. Blue =
  producers consistently *add* to the base at that pixel (positive
  agreement, producers chose noise that lightens the region); red =
  consistently *subtract* (negative agreement, producers chose noise
  that darkens the region).

- **Saturation** encodes `|t|`. Deep color at either end means strong,
  consistent agreement; pale color means weak or inconsistent. The
  colorbar on the right reads in `t` units.

- **`zlim`** is symmetric around zero by default so the neutral color
  aligns with `t = 0`. Pass `zlim = c(-z, z)` to fix the scale across
  panels for direct comparison.

**`palette = "fire"`.** Encodes `|t|` only on a single-hue ramp (pale
yellow at zero -\> deep red at large `|t|`). Use when the question is
*where* producers have a consistent opinion and the direction is not
needed. The `"fire"` view **discards sign by design**; it cannot
distinguish "producers consistently added" from "producers consistently
subtracted". To recover direction at any region of interest, view the
same data with `palette = "diverging"` or pair with
[`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)
of the group-mean CI.

- **Hue intensity** encodes `|t|`. Pale yellow / near-white at low `|t|`
  (so the underlying base face shows through low- agreement regions);
  orange at moderate `|t|`; deep red at large `|t|`. The colorbar reads
  in `|t|` units.

- **`zlim`** defaults to `c(0, max(|t|))` and is asymmetric.

**Common to both palettes.**

- **`threshold`** clips color to the neutral end below
  `|t| < threshold`, making strong-agreement clusters stand out. This is
  descriptive only; it does not provide FWER control. For inferential
  pixel significance, use
  [`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md)
  and render its result directly via `plot(agreement_map_test(...))`, or
  overlay the contours via
  [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md).

- **`base_image`** composites the heatmap on top of a grayscale base
  face so anatomical context shows through. Out-of-mask and subthreshold
  pixels render fully transparent; the per-pixel opacity scales
  `|t| / zlim_max` up to `alpha_max`. Works for both palettes. The color
  bar still shows the full scale so magnitudes are readable off the
  rendered overlay.

- The diverging color convention (blue = positive, red = negative)
  matches
  [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)
  and the cluster-test plots so the same group CI reads consistently
  across the package. The `"fire"` option is unique to this function;
  the CI-overlay and cluster-test plots need to show direction, so they
  do not provide a magnitude-only view.

## See also

[`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)
for the producer-mean counterpart (signed CI, optionally with FWE
contours);
[`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md)
for FWE-controlled significance, and its
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for a
one-call agreement map with contours;
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)
for inferential between-condition tests;
[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
/
[`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
for the optional `mask`.

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

if (FALSE) { # \dontrun{
# Composite the agreement map on the base face for a single
# publication-grade figure. Works for both palettes; the
# "diverging" branch matches plot_ci_overlay()'s color mapping.
sim <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "target",
  signal_region = "eyes", signal_strength = "strong", seed = 1
)
cis <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
plot_agreement_map(cis$signal_matrix,
                   base_image = sim$base_face,
                   threshold  = 2.0,
                   main       = "Agreement t-map over base face")
} # }
```
