# Write rendered CIs to PNG or JPEG files

Renders each column of a `signal_matrix` (one CI) over a base face image
and writes the result to disk. Works for both per-producer matrices (one
column per producer, the `$signal_matrix` returned by
`ci_from_responses_*()`) and group-level matrices (one column per group,
the output of
[`group_ci()`](https://olivethree.github.io/rcisignal/reference/group_ci.md)
or `$group_ci` when the generator is called with `group_by =`).
Filenames are derived from the column names of `signal_matrix`.

Each rendered image uses the same alpha-over compositing path as
[`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)
(default, signed signal with blue=positive, red=negative) or
[`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md)
(`palette = "fire"`, unipolar `|t|`-style yellow-to-red on a per-pixel
CI magnitude).

## Usage

``` r
save_ci_images(
  signal_matrix,
  base_image,
  dir,
  format = c("png", "jpeg"),
  palette = c("diverging", "fire"),
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

  Color palette. `"diverging"` (default; signed CI, blue=positive,
  red=negative, matching
  [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md))
  or `"fire"` (unipolar `|t|`-style, yellow-to-red on `abs(CI)`).

- prefix:

  Optional character scalar overriding the auto-derived filename prefix
  (`"ind_ci_"` for per-producer matrices, `"group_ci_"` for
  group-level). Pass any string to force a custom convention (e.g.
  `prefix = "trust_"`).

- threshold:

  Optional numeric. Pixels with absolute CI value below `threshold` are
  forced to neutral (transparent overlay).

- mask:

  Optional logical vector of length `nrow(signal_matrix)`. Pixels with
  `mask = FALSE` render as base only.

- zlim:

  Optional `c(lo, hi)` color-scale endpoints. When `NULL` (default),
  each image is auto-scaled to its own maximum-absolute value (so colors
  are comparable within a producer / group, not across).

- alpha_max:

  Numeric in `[0, 1]`. Maximum opacity of the heatmap at the color-scale
  top. Default `0.7`.

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

# Per-producer CIs (files: ind_ci_P001.png, ind_ci_P002.png, ...)
save_ci_images(res$signal_matrix, base_image = sim$base_face,
               dir = out)

# Group-level CIs (files: group_ci_A.png, group_ci_B.png)
save_ci_images(res$group_ci, base_image = sim$base_face,
               dir = out, palette = "fire")
} # }
```
