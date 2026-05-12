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
  image's dimensions do not match — useful for catching a
  wrong-resolution mask before it silently corrupts a downstream
  computation.

## Value

Logical vector of length `prod(img_dims)`, column-major.

## See also

[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md),
[`infoval()`](https://olivethree.github.io/rcisignal/reference/infoval.md),
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
