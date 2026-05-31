# Write CIs to PNG or JPEG files

Writes each column of a `signal_matrix` (one CI) to disk as its own
image. The default output matches what
[`rcicr::generateCI()`](https://rdrr.io/pkg/rcicr/man/generateCI.html) /
[`rcicr::generateCI2IFC()`](https://rdrr.io/pkg/rcicr/man/generateCI2IFC.html)
would write for the same CI: a grayscale luminance image of the CI noise
scaled into `[0, 1]` and averaged with the base face (no color palette).
Two palette overlays are available as opt-ins for visualization:
`"diverging"` (matches
[`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md),
signed signal, blue = positive, red = negative) and `"fire"` (matches
`plot_agreement_map(palette = "fire")`, unipolar `|t|`-style
yellow-to-red).

Works for both per-producer matrices (one column per producer, the
`$signal_matrix` returned by `ci_from_responses_*()`) and group-level
matrices (one column per group, the output of
[`group_ci()`](https://olivethree.github.io/rcisignal/reference/group_ci.md)
or `$group_ci` when the generator is called with `group_by =`).
Filenames are derived from the column names of `signal_matrix`.

## Usage

``` r
save_ci_images(
  signal_matrix,
  base_image,
  dir,
  format = c("png", "jpeg"),
  palette = c("grayscale", "diverging", "fire"),
  scaling = c("independent", "constant", "matched", "none"),
  scaling_constant = 0.1,
  prefix = NULL,
  threshold = NULL,
  mask = NULL,
  zlim = NULL,
  alpha_max = 0.7,
  img_dims = NULL,
  quality = 90,
  overwrite = FALSE,
  quiet = FALSE
)
```

## Arguments

- signal_matrix:

  Numeric matrix with non-empty, unique column names. Per-producer or
  group-level; both are accepted.

- base_image:

  Base face image. Either a numeric matrix in `[0, 1]` or a single
  string path to a PNG / JPEG. Used as the underlay for every rendered
  CI.

- dir:

  Output directory. Created (recursively) if missing.

- format:

  Output format. `"png"` (default) or `"jpeg"`.

- palette:

  Color palette. `"grayscale"` (default; raw pixel luminance, matches
  rcicr), `"diverging"` (signed CI on a blue/red ramp, matches
  [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)),
  or `"fire"` (unipolar `|t|`-style yellow-to-red).

- scaling:

  Scaling method for the `"grayscale"` palette, matching rcicr's
  `generateCI(scaling = ...)`: `"independent"` (default; per-CI
  symmetric scaling by `max(|ci|)`), `"constant"` (scale by a
  user-supplied `scaling_constant`, comparable across CIs), `"matched"`
  (range-match each CI to the base image range), or `"none"` (write the
  raw `ci + base` with no scaling; rarely what you want). Ignored when
  `palette != "grayscale"`.

- scaling_constant:

  Numeric constant used when `scaling = "constant"`. Default `0.1`,
  matching rcicr.

- prefix:

  Optional character scalar overriding the auto-derived filename prefix
  (`"ind_ci_"` for per-producer matrices, `"group_ci_"` for
  group-level). Pass any string to force a custom convention (e.g.
  `prefix = "trust_"`).

- threshold:

  Optional numeric. Pixels with absolute CI value below `threshold` are
  forced to 0 (grayscale) or to neutral (palette overlays).

- mask:

  Optional logical vector of length `nrow(signal_matrix)`. Pixels with
  `mask = FALSE` are set to `NA` (grayscale, matching rcicr's
  `applyMask()` semantics) or rendered as base only (palette overlays).

- zlim:

  Optional `c(lo, hi)` color-scale endpoints. Used only for
  `palette = "diverging"` and `palette = "fire"`. Ignored for grayscale.

- alpha_max:

  Numeric in `[0, 1]`. Maximum opacity of the heatmap at the color-scale
  top. Used only for palette overlays. Default `0.7`.

- img_dims:

  Optional integer `c(nrow, ncol)`. Inferred from
  `attr(signal_matrix, "img_dims")` or from a square root of
  `nrow(signal_matrix)`.

- quality:

  JPEG quality in `[0, 100]`. Default `90`. Ignored for PNG.

- overwrite:

  Logical. When `FALSE` (default), the function aborts if any target
  file already exists. When `TRUE`, existing files are silently
  replaced.

- quiet:

  Logical. When `FALSE` (default), emit a one-line `cli` summary at the
  end.

## Value

Invisibly, a character vector of the file paths written.

## Details

Filenames default to `<prefix><colname>.<ext>`, where `<prefix>` is
chosen automatically from `attr(signal_matrix, "ci_level")`:

- `"individual"` (set by `ci_from_responses_*()` on the per-producer
  `$signal_matrix`) -\> `prefix = "ind_ci_"`.

- `"group"` (set by
  [`group_ci()`](https://olivethree.github.io/rcisignal/reference/group_ci.md)
  on its return matrix) -\> `prefix = "group_ci_"`.

- No `ci_level` attribute -\> defaults to `prefix = "ind_ci_"`.

Override the auto-prefix by passing `prefix =` explicitly.

The default `palette = "grayscale"` reproduces rcicr's
`generateCI(..., save_as_png = TRUE)` output exactly: for each CI
column, the raw noise is scaled into `[0, 1]` via the chosen `scaling`
method (default `"independent"`, matching rcicr's default), then
averaged with the base via `(scaled + base) / 2`, then written via
[`png::writePNG()`](https://rdrr.io/pkg/png/man/writePNG.html) (or
[`jpeg::writeJPEG()`](https://rdrr.io/pkg/jpeg/man/writeJPEG.html)) as a
grayscale image. Pass `palette = "diverging"` or `"fire"` instead to
write a colored overlay rendered the same way as the on-screen plot
functions.

## See also

[`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md),
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md),
[`group_ci()`](https://olivethree.github.io/rcisignal/reference/group_ci.md),
[`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md),
[`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md).

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_briefrc_data(
  n_per_condition = 10, n_trials = 60,
  conditions = c("A", "B"), seed = 1
)
res <- ci_from_responses_briefrc(
  sim$data, noise_matrix = sim$noise_matrix,
  base_image = sim$base_face, group_by = "condition"
)

out <- tempfile("ci_export_"); dir.create(out)

# Default: rcicr-style grayscale (raw luminance, no palette).
# Files: ind_ci_P001.png, ind_ci_P002.png, ...
save_ci_images(res$signal_matrix, base_image = sim$base_face,
               dir = out)

# Group-level CIs, same rcicr-style grayscale output.
save_ci_images(res$group_ci, base_image = sim$base_face,
               dir = out)

# Diverging blue/red overlay (rcisignal visualization, not rcicr).
save_ci_images(res$group_ci, base_image = sim$base_face,
               dir = out, palette = "diverging",
               prefix = "diverging_")
} # }
```
