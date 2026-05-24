# Plot an agreement-map test result

Renders the observed per-pixel t-map from
[`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md)
with the same color conventions as
[`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md),
and overlays the FWE-significant pixel boundary as black contours when
`show_contour = TRUE` (the default). Optionally composites the map on a
grayscale base face.

This is the one-call form of the canonical pairing
`plot_agreement_map(signal, ...) + agreement_map_test(...)`: the test
object carries everything the renderer needs (`observed_t`,
`significant_mask`, `img_dims`), so users do not have to re-thread the
source `signal_matrix`.

## Usage

``` r
# S3 method for class 'rcisignal_rel_agreement_map_test'
plot(
  x,
  palette = c("diverging", "fire"),
  threshold = NULL,
  zlim = NULL,
  base_image = NULL,
  alpha_max = 0.7,
  show_contour = TRUE,
  contour_col = "black",
  contour_lwd = 1,
  main = "Agreement t-map (FWE contours)",
  ...
)
```

## Arguments

- x:

  A
  [`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md)
  result.

- palette:

  `"diverging"` (default; signed t, blue = positive, red = negative) or
  `"fire"` (`|t|` on a single-hue ramp). See
  [`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md)
  for the full Reading-the-plot discussion; the same conventions apply
  here.

- threshold:

  Optional positive numeric. Pixels with `|t| < threshold` render as the
  neutral color (descriptive only; FWE control is already in
  `significant_mask`).

- zlim:

  Numeric `c(low, high)` for the color scale. Defaults to a symmetric
  `c(-max|t|, max|t|)` for diverging or `c(0, max|t|)` for fire.

- base_image:

  Optional. Numeric matrix or path to PNG/JPEG. When supplied, the t-map
  is composited on top of the grayscale base; out-of-mask and
  subthreshold pixels render fully transparent.

- alpha_max:

  Numeric in `[0, 1]`. Maximum opacity at the color-scale top when
  `base_image` is supplied. Default 0.7.

- show_contour:

  Logical. Draw the FWE-significant pixel boundary as black contours on
  top of the t-map. Default `TRUE`.

- contour_col, contour_lwd:

  Significance-contour color and line width.

- main:

  Plot title.

- ...:

  Reserved for future use.

## Value

Invisibly the input `x`.

## See also

[`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md),
[`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md),
[`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md).
