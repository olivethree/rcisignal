# Read an image-based face-region mask from a PNG or JPEG

Reads a binary mask image from disk and returns a logical vector of
length `prod(img_dims)` in column-major order, the format the rest of
the package expects. White / light pixels (above `threshold`) become
`TRUE`, dark pixels become `FALSE`. Use this for image-based masks
created with
[webmorphR::mask_oval()](https://debruine.github.io/webmorphR/), painted
in GIMP, drawn with PowerPoint shapes, or any other tool that produces a
binary PNG/JPEG.

Companion to
[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
(parametric oval). Either function returns the same logical-vector
shape.

## Usage

``` r
read_face_mask(path, threshold = 0.5, invert = FALSE, expected_dims = NULL)
```

## Arguments

- path:

  Path to a PNG or JPEG mask image.

- threshold:

  Numeric in `[0, 1]`. Pixels with luminance above this become `TRUE`.
  Default `0.5` (mid-grey).

- invert:

  If `TRUE`, the mask is inverted (useful for black-on-white masks).
  Default `FALSE`.

- expected_dims:

  Optional integer `c(nrow, ncol)`. When set, aborts if the loaded
  image's dimensions do not match. Useful for catching a
  wrong-resolution mask before it silently corrupts a downstream
  computation.

## Value

Logical vector of length `prod(img_dims)`, column-major.

## See also

[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
for the parametric oval / sub-region constructor. Mask consumers:
[`infoval()`](https://olivethree.github.io/rcisignal/reference/infoval.md),
[`pixel_t_test()`](https://olivethree.github.io/rcisignal/reference/pixel_t_test.md),
[`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md),
[`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md),
[`rel_icc()`](https://olivethree.github.io/rcisignal/reference/rel_icc.md),
[`rel_split_half()`](https://olivethree.github.io/rcisignal/reference/rel_split_half.md),
[`rel_loo()`](https://olivethree.github.io/rcisignal/reference/rel_loo.md),
[`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md),
[`run_reliability()`](https://olivethree.github.io/rcisignal/reference/run_reliability.md),
[`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md),
[`run_discriminability_pairwise()`](https://olivethree.github.io/rcisignal/reference/run_discriminability_pairwise.md),
[`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md),
[`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md).
Plot helpers:
[`plot_face_mask()`](https://olivethree.github.io/rcisignal/reference/plot_face_mask.md),
[`plot_mask_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_mask_overlay.md).

## Examples

``` r
if (FALSE) { # \dontrun{
fm <- read_face_mask("masks/oval_256.png",
                     expected_dims = c(256L, 256L))
iv <- infoval(signal_matrix, noise_matrix, trial_counts,
              mask = fm, iter = 1000L, seed = 1L)
} # }
```
