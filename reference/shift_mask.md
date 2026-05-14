# Shift a face-region mask by a number of pixels

Slides a logical mask by `down` and `right` pixels and returns the
shifted mask in the same shape as the input. Convenience wrapper around
the manual recipe shown in vignette section 4.5 for tuning an
*elliptical*
[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
region (`"full"`, `"nose"`, `"mouth"`, `"upper_face"`, `"lower_face"`)
on a base image whose features sit a few pixels off the default
geometry.

For the *rectangle* regions (`"eyes"`, `"left_eye"`, `"right_eye"`),
prefer the `region_bounds` argument of
[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md):
passing `region_bounds = c(row_min, row_max, col_min, col_max)` in 0-1
image fractions is more precise than counting pixels and is the
canonical tuning route for those regions.

Pixels shifted off the image are dropped, and the newly exposed edge is
filled with `FALSE` (i.e., outside the mask).

## Usage

``` r
shift_mask(mask, down = 0L, right = 0L, img_dims = NULL)
```

## Arguments

- mask:

  Logical vector of length `prod(img_dims)` (column-major, as returned
  by
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md))
  or a logical matrix.

- down, right:

  Integer pixel offsets. Positive `down` moves the mask towards the
  bottom of the image; positive `right` moves it towards the right edge.
  Negative values move up / left. Defaults are `0` (no shift).

- img_dims:

  Integer `c(nrow, ncol)` (or a single integer for a square image).
  Required when `mask` is a vector; ignored when `mask` is already a
  matrix.

## Value

Logical mask of the same shape as the input (matrix-in, matrix-out;
vector-in, vector-out).

## See also

[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
(the `region_bounds` route for the rectangle regions),
[`plot_mask_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_mask_overlay.md),
[`plot_face_mask()`](https://olivethree.github.io/rcisignal/reference/plot_face_mask.md).

## Examples

``` r
# Reshape the default mouth mask, slide it down 20 pixels on a
# 256-pixel image (about 8% of the height), and confirm it
# still has the right number of TRUE pixels.
m <- make_face_mask(c(256L, 256L), region = "mouth")
shifted <- shift_mask(m, down = 20L, img_dims = c(256L, 256L))
identical(sum(m), sum(shifted))
#> [1] TRUE

# Matrix-in, matrix-out.
mat <- matrix(m, 256L, 256L)
shifted_mat <- shift_mask(mat, down = 20L, right = 8L)
dim(shifted_mat)
#> [1] 256 256
```
