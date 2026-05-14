# Convert pixel-coordinate rectangle bounds to 0-1 image fractions

Bridge between visual inspection of a base image (where people naturally
think in pixel rows and columns) and the `region_bounds` argument of
[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md),
which expects 0-1 image fractions. Returns a length-4 numeric vector
ready to pass straight to `make_face_mask(region_bounds = ...)`.

## Usage

``` r
region_bounds_from_pixels(row_min, row_max, col_min, col_max, img_dims)
```

## Arguments

- row_min, row_max, col_min, col_max:

  Integer pixel coordinates of the rectangle edges (1-indexed, with row
  1 at the top of the image and col 1 at the left). Must satisfy
  `row_min <= row_max` and `col_min <= col_max`, and all four values
  must lie inside the image.

- img_dims:

  Integer `c(nrow, ncol)` of the target image (or a single integer for a
  square image).

## Value

Numeric length-4 vector `c(row_min, row_max, col_min, col_max)` in 0-1
image fractions.

## See also

[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md).

## Examples

``` r
# On a 256x256 base, "rows 110-140, cols 60-115" becomes:
region_bounds_from_pixels(
  row_min = 110, row_max = 140,
  col_min = 60,  col_max = 115,
  img_dims = c(256L, 256L)
)
#> [1] 0.4274510 0.5450980 0.2313725 0.4470588

# Pass straight to make_face_mask():
bounds <- region_bounds_from_pixels(110, 140, 60, 115,
                                    img_dims = 256L)
m <- make_face_mask(c(256L, 256L), region = "left_eye",
                    region_bounds = bounds)
mean(m)
#> [1] 0.02648926
```
