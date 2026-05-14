# Build an oval face-region mask for a square image

Returns a logical vector of length `prod(img_dims)` marking a face
region (or a face sub-region) centred on the image. Pass to
[`infoval()`](https://olivethree.github.io/rcisignal/reference/infoval.md)
to restrict both observed and reference Frobenius norms to the masked
region; pass to any `rel_*()` function via row-subsetting
(`signal_matrix[mask, ]`) to compute reliability or discriminability on
a single anatomical region.

Eight regions are supported:

- `"full"` (default): the full face oval. Defaults follow the geometry
  used by Schmitz, Rougier, & Yzerbyt (2024) and tunable via `centre` /
  `half_width` / `half_height`.

- `"eyes"`: a wide rectangle covering both eyes ear-to-ear and from the
  eyebrows down to just below the eye line. Independent of the full-oval
  geometry; tune via `region_bounds`.

- `"left_eye"`, `"right_eye"`: smaller rectangles around the viewer's
  left and right eye respectively. Independent of the full-oval
  geometry; tune via `region_bounds`.

- `"nose"`: a narrow vertical ellipse along the midline.

- `"mouth"`: a wide-and-short ellipse below centre.

- `"upper_face"`, `"lower_face"`: top and bottom halves of the full face
  oval.

Sub-region geometries are heuristic approximations matched to a typical
centred face on a square base (e.g. 256x256). For non-default base
images, the elliptical regions (`"nose"`, `"mouth"`) scale with `centre`
/ `half_width` / `half_height`; the three rectangle regions (`"eyes"`,
`"left_eye"`, `"right_eye"`) ignore those and are tuned directly via
`region_bounds`.

## Usage

``` r
make_face_mask(
  img_dims,
  region = c("full", "eyes", "left_eye", "right_eye", "nose", "mouth", "upper_face",
    "lower_face"),
  centre = c(0.5, 0.5),
  half_width = 0.35,
  half_height = 0.45,
  region_bounds = NULL
)
```

## Arguments

- img_dims:

  Integer `c(nrow, ncol)`.

- region:

  Character. One of `"full"`, `"eyes"`, `"left_eye"`, `"right_eye"`,
  `"nose"`, `"mouth"`, `"upper_face"`, `"lower_face"`.

- centre:

  Numeric `c(row, col)` in the 0-1 coordinate range. Default
  `c(0.5, 0.5)`. Used by the elliptical regions only; the rectangle eye
  regions ignore this argument.

- half_width:

  Full-face ellipse horizontal half-axis, as a fraction of image width.
  Default 0.35. Used by the elliptical regions only.

- half_height:

  Full-face ellipse vertical half-axis, as a fraction of image height.
  Default 0.45. Used by the elliptical regions only.

- region_bounds:

  Numeric vector `c(row_min, row_max, col_min, col_max)` in the 0-1
  image fraction range, used to override the default rectangle bounds
  for the three rectangle regions (`"eyes"`, `"left_eye"`,
  `"right_eye"`). Each pair must be ordered (`row_min < row_max`,
  `col_min < col_max`) and lie within `[0, 1]`. Default `NULL` uses the
  built-in defaults. Errors when supplied for a non-rectangle region.

## Value

Logical vector of length `prod(img_dims)`, column-major to match the
package's image vectorisation convention.

## References

Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the brief
reverse correlation: an improved tool to assess visual representations.
*European Journal of Social Psychology*.
[doi:10.1002/ejsp.3100](https://doi.org/10.1002/ejsp.3100)

## See also

[`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md),
[`infoval()`](https://olivethree.github.io/rcisignal/reference/infoval.md),
[`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md),
[`plot_face_mask()`](https://olivethree.github.io/rcisignal/reference/plot_face_mask.md),
[`plot_mask_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_mask_overlay.md)

## Examples

``` r
full      <- make_face_mask(c(128L, 128L))
eyes      <- make_face_mask(c(128L, 128L), region = "eyes")
left_eye  <- make_face_mask(c(128L, 128L), region = "left_eye")
right_eye <- make_face_mask(c(128L, 128L), region = "right_eye")
mouth     <- make_face_mask(c(128L, 128L), region = "mouth")
c(full = mean(full), eyes = mean(eyes),
  left_eye = mean(left_eye), right_eye = mean(right_eye),
  mouth = mean(mouth))
#>       full       eyes   left_eye  right_eye      mouth 
#> 0.48754883 0.17395020 0.02062988 0.02062988 0.01440430 

# Tune a rectangle region by passing explicit bounds in 0-1
# image fractions: c(row_min, row_max, col_min, col_max).
left_eye_tuned <- make_face_mask(
  c(256L, 256L), region = "left_eye",
  region_bounds = c(0.40, 0.50, 0.24, 0.44)
)
```
