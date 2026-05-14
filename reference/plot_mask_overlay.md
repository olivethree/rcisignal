# Sanity-check a face mask against a base image

Renders a base face with a translucent mask overlay so you can verify
the mask aligns with the anatomy you intended before feeding it into
[`infoval()`](https://olivethree.github.io/rcisignal/reference/infoval.md),
[`run_reliability()`](https://olivethree.github.io/rcisignal/reference/run_reliability.md),
or any other downstream analysis. Avoids the usual hand-rolled
`image() + adjustcolor()` workflow.

Companion to
[`plot_face_mask()`](https://olivethree.github.io/rcisignal/reference/plot_face_mask.md)
(mask alone, optional base) and
[`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)
(signal + base): use `plot_mask_overlay()` when the question is "does
the mask cover the right region of this specific base image?".

## Usage

``` r
plot_mask_overlay(
  base_image,
  mask = NULL,
  region = NULL,
  region_bounds = NULL,
  alpha = 0.35,
  overlay_col = "red",
  main = NULL
)
```

## Arguments

- base_image:

  Either a numeric matrix in `[0, 1]` (grayscale, `nrow x ncol`) or a
  single character path to a PNG/JPEG file. PNG/JPEG paths are read via
  the same grayscale conversion the rest of the package uses (BT.709
  luminance for RGB inputs).

- mask:

  One of: a logical or numeric vector of length `prod(dim(base_image))`
  (column-major, matching the convention used by
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  /
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md));
  a logical or numeric matrix of the same dimensions as `base_image`; or
  a single character path to a PNG/JPEG mask file (resolved via
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)).
  Numeric inputs are thresholded at `0.5`. Pass `NULL` (the default) and
  a `region` to build the mask internally via
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md).

- region:

  Optional character region name for the convenience shortcut: passing
  e.g. `region = "left_eye"` builds the mask via
  `make_face_mask(dim(base_image), region, region_bounds)` and skips the
  separate construction step. Mutually exclusive with `mask`.

- region_bounds:

  Optional length-4 numeric vector forwarded to
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  when `region` is one of the rectangle regions (`"eyes"`, `"left_eye"`,
  `"right_eye"`). Ignored otherwise.

- alpha:

  Numeric in `[0, 1]`. Opacity of the mask overlay. Default `0.35`.

- overlay_col:

  Single colour for the mask overlay. Default `"red"`.

- main:

  Optional plot title (character of length 1). Default `NULL` (no
  title).

## Value

Invisibly `NULL`. Called for the side-effect of drawing on the active
graphics device.

## See also

[`plot_face_mask()`](https://olivethree.github.io/rcisignal/reference/plot_face_mask.md),
[`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md),
[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md),
[`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Parametric oval mask plus a flat synthetic base.
n_side <- 128L
base   <- matrix(0.5, n_side, n_side)
eyes   <- make_face_mask(c(n_side, n_side), region = "eyes")
plot_mask_overlay(base, eyes)
} # }

if (FALSE) { # \dontrun{
# Use the simulator's base face and a sub-region to verify
# that the parametric mask lands on the intended anatomy.
sim   <- simulate_briefrc_data(
  n_per_condition = 5, n_trials = 10, conditions = "target",
  seed = 1
)
mouth <- make_face_mask(dim(sim$base_face), region = "mouth")
plot_mask_overlay(sim$base_face, mouth, main = "Mouth region")
} # }

if (FALSE) { # \dontrun{
# Mask loaded from a PNG file; pass both as paths.
plot_mask_overlay("path/to/base.png", "path/to/mask.png")
} # }

# Region shortcut: skip the make_face_mask() call.
base <- matrix(0.5, 128L, 128L)
plot_mask_overlay(base, region = "left_eye")
```
