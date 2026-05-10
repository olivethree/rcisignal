# Render a group CI as a translucent agreement-map overlay on a base image

Plots the producer-mean classification image (or the per-pixel *t*-map
from
[`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md))
as a diverging-palette translucent layer composited over the base face
image. Useful as a publication-grade single figure that combines
anatomical context (the base) with the empirical signal (the CI /
agreement map). Optionally overlays significance contours from an
[`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md)
result.

Pure base graphics — no `grid` / `ggplot2` dependency. Returns invisibly
the composed `nrow x ncol x 3` numeric raster so it can be saved to disk
via [`png::writePNG()`](https://rdrr.io/pkg/png/man/writePNG.html) /
[`jpeg::writeJPEG()`](https://rdrr.io/pkg/jpeg/man/writeJPEG.html) if
needed.

## Usage

``` r
plot_ci_overlay(
  signal_matrix,
  base_image,
  img_dims = NULL,
  test = NULL,
  mask = NULL,
  threshold = NULL,
  alpha_max = 0.7,
  palette = c("diverging"),
  contour_col = "black",
  contour_lwd = 1,
  main = NULL
)
```

## Arguments

- signal_matrix:

  Pixels x participants; the producer mean is computed and rendered as
  the heatmap layer. Alternatively, a numeric vector of length
  `prod(img_dims)` (e.g. `agreement_map_test(...)$observed_t`) which is
  rendered directly.

- base_image:

  Either a numeric matrix (`nrow x ncol`, grayscale, values in 0-1) or a
  path to a PNG/JPEG file. When it's a path, the image is read as
  grayscale (BT.709 luminance for RGB inputs) the same way as
  [`read_cis()`](https://olivethree.github.io/rcisignal/reference/read_cis.md).

- img_dims:

  Integer `c(nrow, ncol)`. If `NULL`, inferred from
  `attr(signal_matrix, "img_dims")` or, when `signal_matrix` is a
  vector, the dimensions of `base_image`.

- test:

  Optional
  [`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md)
  result. When supplied, significance contours (around
  `significant_mask`) are drawn on top of the heatmap.

- mask:

  Optional logical vector of length `prod(img_dims)` restricting the
  visible overlay region. Pixels outside the mask render as the bare
  base image. Build with
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  (parametric oval and sub-regions) or
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
  (PNG/JPEG mask file). The mask is column-major to match the package's
  image vectorisation convention.

- threshold:

  Optional numeric. Pixels with `|signal| < threshold` are rendered as
  fully transparent (only the base shows through). Default `NULL` (no
  threshold).

- alpha_max:

  Numeric in `[0, 1]`. Maximum opacity at the peak of the signal
  magnitude. Default 0.7.

- palette:

  Character. Currently only `"diverging"` is implemented
  (red-white-blue). Future-proofed by argument.

- contour_col, contour_lwd:

  Significance-contour colour and line width when `test` is supplied.
  Defaults: `"black"`, `1.0`.

- main:

  Optional plot title.

## Value

Invisibly the composed `nrow x ncol x 3` raster. The plot is drawn on
the active device as a side effect.

## Reading the plot

- **Colour** encodes the sign of the producer-mean signal at each pixel.
  Blue = positive (producers' average mask is brighter than the base at
  that pixel); red = negative (darker than the base); pixels rendered as
  the bare base = at or near zero.

- **Opacity** encodes the magnitude of the signal, scaled to `alpha_max`
  at the global peak `|signal|`. Faint colour means weak agreement;
  saturated colour means a strong, consistent producer-mean deflection.

- **Black contours** (only drawn when `test` is supplied) trace the
  boundary of the FWE-significant pixel set returned by
  [`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md)
  at its `alpha`. Pixels inside the contour are individually significant
  under the max-\|t\| null; pixels outside are not.

- Colour convention matches
  [`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md)
  and the cluster-test plots so the same group CI reads consistently
  across the package.

## See also

[`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md),
[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal call-signature demo with synthetic inputs and a flat base.
n_side <- 32L
n_pix  <- n_side * n_side
n_prod <- 20L
set.seed(1)
signal_matrix <- matrix(rnorm(n_pix * n_prod), n_pix, n_prod)
base_image    <- matrix(0.5, n_side, n_side)
plot_ci_overlay(signal_matrix, base_image,
                img_dims = c(n_side, n_side))
} # }

if (FALSE) { # \dontrun{
# Same function, richer input: plant a signal in the eye region and
# overlay it on the simulator's base face. The eye region should
# light up against the anatomical context.
sim <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "target",
  signal_region = "eyes", signal_strength = "strong", seed = 1
)
cis  <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
plot_ci_overlay(cis$signal_matrix, sim$base_face)
} # }

if (FALSE) { # \dontrun{
# Same overlay with FWE-significant pixels outlined: feed an
# agreement_map_test() result via `test = `. Black contours trace
# the boundary of pixels significant at the test's alpha.
sim   <- simulate_briefrc_data(
  n_per_condition = 20, n_trials = 60, conditions = "target",
  signal_region = "eyes", signal_strength = "strong", seed = 1
)
cis   <- ci_from_responses_briefrc(sim$data, noise_matrix = sim$noise_matrix)
agree <- agreement_map_test(cis$signal_matrix,
                            n_permutations = 500L, seed = 1)
plot_ci_overlay(cis$signal_matrix, sim$base_face,
                test = agree,
                main = "CI overlay + significant pixels")
} # }
```
