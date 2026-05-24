# rcisignal 0.2.0

## Breaking changes

* `group_ci()` redesigned. The vector / list-of-vectors `by`
  argument is gone. New signature:

  ```r
  group_ci(signal_matrix, responses, by,
           col_participant = "participant_id", drop = TRUE)
  ```

  `responses` is the same trial-level data frame you handed to
  `ci_from_responses_*()`; `by` is the name of the grouping
  column (or a character vector of column names for factorial
  grouping). Producer-to-group alignment happens internally via
  `colnames(signal_matrix)`. Migration:

  ```r
  # Before:
  cond <- responses$condition[match(colnames(sm), responses$participant_id)]
  gcis <- group_ci(sm, by = cond)

  # Now:
  gcis <- group_ci(sm, responses, by = "condition")
  ```

  `group_ci()` also validates that every producer in
  `colnames(signal_matrix)` is present in `responses`, and that
  each producer's `by` column(s) are constant across their
  rows. Both checks fail with teaching messages naming the
  offending producer.

* `shift_mask()` arguments renamed and `vertical` sign flipped.
  `down` -> `vertical`, `right` -> `horizontal`. The new
  `vertical` follows the math / y-axis-up idiom (positive moves
  the mask up, negative moves it down); `horizontal` keeps the
  intuitive convention (positive right, negative left). Migration:
  rename `down = N` -> `vertical = -N` (note the sign flip),
  `right = N` -> `horizontal = N`. No deprecation alias.

* Column-name arguments standardized on the `col_*` convention
  in the two CI builders. Renamed:
  `ci_from_responses_briefrc(participant_col, stimulus_col, response_col)`
  -> `(col_participant, col_stimulus, col_response)`; same for
  `ci_from_responses_2ifc()`. Every other responses-consuming
  function (`check_*`, `diagnose_*`, `compute_*`, `run_*`)
  already used `col_*`; the package is now uniform. No
  deprecation alias; calls using the old names fail with the
  standard R unused-argument error.

## New features

* `infoval()` gains a `responses = NULL`, `col_participant =
  "participant_id"` path. When `trial_counts` is omitted,
  `infoval()` derives it internally via
  `table(responses[[col_participant]])` and matches against
  `colnames(signal_matrix)`. Existing scripts that pass
  `trial_counts` directly continue to work unchanged.

* `run_reliability()` validates `noise_matrix` at function
  entry when `null = "random_responders"`, with a teaching
  message pointing at `read_noise_matrix()`. Replaces a deeper
  cryptic abort that used to fire inside `rel_split_half()`.

## Documentation

* Vignette §1.3 (two-stage pattern), §10.5, and the
  worked-example cross-references updated to use the new
  `group_ci()` idiom. §9 (infoval) and roxygen `@examples`
  on `infoval()` updated to lead with the `responses` path.

# rcisignal 0.1.9

## Breaking change

* `ci_from_responses_2ifc()` argument `baseimage` renamed to
  `base_image` to match `ci_from_responses_briefrc()`. The new
  argument now accepts the same three input forms as Brief-RC:
  a numeric matrix in `[0, 1]`, a path to a PNG / JPEG, or a
  string label naming an entry in the rdata's `base_faces` list
  (the historical 2IFC form). Matrix and path inputs are injected
  into a fresh temporary copy of the rdata under a synthetic
  label before the call into rcicr. No deprecation alias; calls
  using `baseimage = "..."` will fail with the standard R
  unused-argument error.

* The same rename cascades through the five downstream 2IFC
  diagnostic functions and the noise reader:
  `diagnose_infoval()`, `compute_infoval_summary()`,
  `check_rt_infoval_consistency()`, `check_response_inversion()`,
  `run_diagnostics()`, plus `read_noise_matrix(base_image = ...)`.
  All accept the same three forms via the same internal resolver.

## Documentation

* Tutorial vignette cleaned up for readability per a user-driven
  audit:
  - All user-facing prose (vignette, README, NEWS, roxygen, `cli`
    message bodies) switched from British to American English.
    Argument names and `match.arg()` value strings (e.g.,
    `centre`, `colour_bar`, `method = "normalised"`) were left
    in place to preserve the documented API surface.
  - `DESCRIPTION` `Language` field updated `en-GB` -> `en-US`.
  - `plot_ci_correlogram()`, `plot_ci_distance_matrix()`, and
    `plot_ci_mds()` `@param cis` and the worked-example chunks
    in vignette §10.5 / §10.6 / §10.7 now lead with the named-
    matrix input form (one column per CI, built outside the call
    with `cbind(name = rowMeans(...))`). The named-list and
    `group_ci()` forms are still accepted; they are documented
    as the back-compatibility / convenience paths. Mixed
    group-CI + individual-producer correlogram example added.
  - §6.1 (compute CIs) now cross-links to §1.3 (two-stage
    pattern) and `group_ci()` as the canonical group-CI route.
  - Removed developer-facing asides that leaked into the
    user-guide register (the "eval = FALSE for faster vignette
    rendering" note in §12.0; the precompute-null reuse aside
    in §7.1; the `attr(res$signal_matrix, "source")` and
    `attr(signal, "source")` example lines in §6.1 and §6.2).

# rcisignal 0.1.8

## New features

* `group_ci()` (new exported function). Stage-2 aggregator that
  collapses a per-producer `signal_matrix` (pixels x n_producers,
  as returned by `ci_from_responses_briefrc()` /
  `ci_from_responses_2ifc()`) into a per-group matrix (pixels x
  n_groups) suitable for the existing distance-matrix, MDS, and
  correlogram plot functions. Accepts an atomic grouping vector
  or a named list of vectors (factorial grouping). Returns a
  numeric matrix classed `c("rcisignal_group_ci", "matrix",
  "array")` with per-group producer counts in `attr(., "n")` and
  `img_dims` inherited from the input. Replaces the hand-rolled
  `rowMeans(signal_matrix[, group_idx])` users were writing.

* New S3 class `rcisignal_group_ci` with `print()` and
  `as.list()` methods. The print method surfaces per-group
  producer counts and reminds the reader that per-producer
  information has been averaged out.

## Robustness

* The stage-1 functions now error with a teaching message when
  handed a `group_ci()` result instead of a per-producer
  `signal_matrix`. The message names the two-stage pattern and
  points at `vignette("rcisignal")`. Affected functions:
  `infoval()`, `rel_split_half()`, `rel_icc()`, `rel_loo()`,
  `rel_cluster_test()` (both `_a` / `_b`), `rel_dissimilarity()`
  (both), `agreement_map_test()`, `pixel_t_test()` (both),
  `run_reliability()`, `run_discriminability()` (both), and
  `run_discriminability_pairwise()` (any list element). The
  internal `group_mean_z()` gets the same guard. Pre-existing
  behavior for legitimate per-producer matrices is unchanged.

## Documentation

* New "Two-stage pattern" section in `vignette("rcisignal")`,
  placed immediately after the Overview. The package now names
  the two stages: stage 1 (per-producer `signal_matrix`,
  produced by `ci_from_responses_*()`, consumed by every
  reliability / discriminability / infoVal function) and stage 2
  (group-averaged matrix, produced by `group_ci()`, consumed by
  the distance-matrix / MDS / correlogram plots).

* pkgdown reference index regrouped to mirror the two stages.
  Visitors to the package website now see "Stage 1 —
  per-producer CIs" and "Stage 2, group CIs" as the first two
  reference sections.

* README gains a small two-stage diagram showing how `responses`
  + `noise_matrix` + `base_image` flow through stage 1 into the
  reliability / infoVal / discriminability surfaces, and
  optionally into stage 2 via `group_ci()` for RDM / MDS /
  correlogram analyses.

* `@seealso` on `ci_from_responses_briefrc()` and
  `ci_from_responses_2ifc()` now points at `group_ci()` as the
  optional stage-2 follow-on.

# rcisignal 0.1.7

## New features

* `plot_agreement_map()` gains `base_image` and `alpha_max`
  arguments. When `base_image` is supplied (numeric matrix or
  PNG/JPEG path), the t-map is composited on top of the grayscale
  base; out-of-mask and subthreshold pixels render fully
  transparent, and the per-pixel opacity scales `|t| / zlim_max`
  up to `alpha_max`. Works for both `palette = "diverging"` and
  `palette = "fire"`. With `base_image = NULL` (default) the
  existing flat-panel rendering is unchanged.

* `plot.rcisignal_rel_cluster_test()` gains the same `base_image`
  and `alpha_max` arguments, so the between-condition t-map (or
  TFCE map) can be rendered with anatomical context in one call.
  Diverging palette only; the cluster t-map is intrinsically
  signed and is not given a `|t|`-only view.

* `plot.rcisignal_rel_pairwise_report()` gains `base_image` and
  `alpha_max`, threaded through to each per-pair cluster-test
  plot. When `base_image` is a path string it is read once and
  reused for every panel.

* `plot.rcisignal_rel_agreement_map_test()` (new S3 method).
  `plot()` of an `agreement_map_test()` result now renders the
  observed t-map directly, with the FWE-significant pixel
  boundary outlined in black by default. Accepts `palette`,
  `threshold`, `zlim`, `base_image`, `alpha_max`, and a
  `show_contour` toggle. To support this, the
  `agreement_map_test()` result gains an `$img_dims` field
  (additive; existing fields unchanged).

## Documentation

* Reciprocal `@seealso` cross-links between `plot_agreement_map()`
  and `plot_ci_overlay()`; both now point at each other and at
  `agreement_map_test()`. Vignette §10 reframed around the three
  composable surfaces (`plot_agreement_map()` /
  `plot_ci_overlay()` / cluster-test plot) so the symmetric
  `base_image` workflow is visible from the start of the section.

* `make_face_mask()` and `read_face_mask()` `@seealso` blocks
  expanded to list every mask-accepting consumer
  (`infoval()`, `pixel_t_test()`, `agreement_map_test()`,
  `rel_cluster_test()`, `rel_icc()`, `rel_split_half()`,
  `rel_loo()`, `rel_dissimilarity()`, `run_reliability()`,
  `run_discriminability()`, `run_discriminability_pairwise()`,
  `plot_agreement_map()`, `plot_ci_overlay()`). Discoverability
  is now complete from either constructor.

## Internal

* Private helpers `resolve_base_for_overlay()` and
  `composite_rgb_over_gray()` factored into `R/utils.R` so
  `plot_ci_overlay()`, `plot_agreement_map()`, and the
  cluster-test plot method share one base-loading and
  compositing path. No user-visible signature change.

# rcisignal 0.1.6

## New features

* `plot_agreement_map()` gains a unipolar `palette = "fire"` option
  that displays `|t|` on a single-hue (pale yellow to deep red)
  ramp. Use when the question is "where do producers have a
  consistent opinion" and direction is not needed; the `"fire"`
  view discards sign by design. The default `palette = "diverging"`
  is unchanged in API; see the color-direction fix below.
* `plot_ci_correlogram()` (new exported function). Renders a
  publication-ready Pearson-`r` matrix across multiple group-mean
  CIs. Inputs: a named list of CIs (vectors or per-producer
  signal matrices; group means computed automatically). Options:
  full / upper / lower triangle, optional face / upper-face /
  lower-face masking via `make_face_mask()`, three diverging
  palettes (RdBu / PuOr / BrBG), direct save to PNG (600 dpi) or
  PDF via the `file = ` argument with extension auto-detection.
  The §12.6 worked-example figure is now produced by this
  function.

* `plot_ci_distance_matrix()` (new exported function). All-vs-all
  pairwise Euclidean distance matrix across a named list of
  group-mean CIs. Same beginner-friendly input format as
  `plot_ci_correlogram()`; the magnitude metric recommended by
  §8.3 instead of Pearson `r`. Supports raw (`method = "raw"`,
  default) or resolution-normalised (`method = "normalised"`)
  distance, optional face / upper-face / lower-face masking, four
  sequential colorblind-safe palettes (viridis / inferno /
  plasma / rocket), full / upper / lower triangle rendering, and
  direct PNG (600 dpi) / PDF save.

* `plot_ci_mds()` (new exported function). Classical
  multidimensional-scaling projection of a named list of
  group-mean CIs (`stats::cmdscale()`). Auto-selects the smallest
  dimensionality whose Kruskal stress-1 against the original
  Euclidean distances reaches the "good" band (default threshold
  `0.05`, Kruskal 1964); renders a grid of all `choose(k, 2)`
  pairwise dimension panels when the auto-selected `k > 2`.
  Users can force a specific `k` (e.g., `k = 2L` for a single
  paper-figure panel). Returns an S3 object (`rcisignal_mds`)
  with a one-screen `print()` method; `$mds_points` exposes the
  `n_cis x k_selected` coordinate matrix in the Euclidean MDS
  space, `$stress_by_k` and `$variance_pct_by_k` expose the
  full dimensionality-selection trace. Supports optional
  categorical `groups` (point color) and `shapes` (point pch)
  arguments for grouped scatters in multi-condition designs.
  Same masking and PNG/PDF save options as the other plot
  helpers.

* New private helper `prepare_ci_matrix()` (internal, not
  exported) centralises the named-list-of-CIs validation and
  per-producer reduction across `plot_ci_correlogram()`,
  `plot_ci_distance_matrix()`, and `plot_ci_mds()`. No
  user-visible behavior change to `plot_ci_correlogram()`.

## Behavioral change

* Color direction fixed on the diverging-palette plots that use
  `hcl.colors("RdBu", ...)`. `plot_agreement_map(palette =
  "diverging")` and `plot.rcisignal_rel_cluster_test()` previously
  rendered positive values as red and negative as blue, opposite
  to what their help pages and to `plot_ci_overlay()` documented.
  Both now render positive = blue, negative = red, matching every
  other diverging plot in the package and matching the
  long-standing help-page docs. Visual effect: rendered figures
  from prior versions of these two functions look like their
  colors have been swapped. The new `plot_ci_correlogram()`
  follows the same convention from the start.
* The worked-example correlogram figure
  (`vignettes/figures/oliveira_2019/trait_ci_correlogram.png`)
  is regenerated to match the corrected convention. Same data,
  flipped colors.

## Breaking change

* `plot_agreement_map(palette = "viridis")` is removed. The
  previous viridis branch plotted signed `t` against a palette
  with no neutral midpoint, so zero rendered at a perceptually
  arbitrary color. The new `"fire"` option supersedes that use
  case correctly. The argument is now
  `palette = c("diverging", "fire")`.

## Documentation

* `plot_agreement_map()` help page rewritten to be explicit that
  on the default diverging palette both deep red and deep blue
  indicate **strong** agreement (sign in hue, magnitude in
  saturation); "no agreement" is the neutral color (white), not
  red. The new `"fire"` palette is documented as a sign-discarding
  magnitude-only view with explicit pointers to recover direction
  via the diverging palette or `plot_ci_overlay()`.
* Vignette §10.2 gains a one-paragraph aside introducing the
  fire option and reiterating the two-channel reading of the
  diverging palette.
* `rel_dissimilarity()`: the Pearson correlation fields
  (`$correlation`, `$boot_cor`, `$ci_cor`, `$boot_se_cor`) are
  reframed as a **secondary** summary instead of being slated
  for removal. They stay in the API. The help page, `print()`
  output, plot panel title ("Pearson r (secondary)"), and
  vignette §8.3 now explain (a) why Euclidean distance is the
  recommended primary statistic (image-domain scaffolding gives
  `r` a positive chance baseline that does not cleanly mean
  "similar"), and (b) how to use `r` carefully if it must be
  reported (relative comparisons across pairs against a
  permutation null, not absolute values against zero). The prior
  "will be removed in v0.2.0" wording is withdrawn.
* Vignette: editorial sweep of dev-facing prose. Replaced
  "Loaded from cache, ... on this dataset:" lead-ins (five
  occurrences) with neutral "On this dataset, ... :" framing.
  Removed roadmap aside about a possible future Rcpp accelerator
  for `simulate_*_data()`. Condensed the §4.2 list of rdata
  bookkeeping fields. Replaced two "self-critical note" author
  asides in §8.2 and §8.3 with neutral "one caveat" framing.
  Dropped a `v1.0.x` version label in §12.2. Rephrased the §14
  pointer at NEWS.md as `news(package = "rcisignal")`.

## Internal

* `plot_agreement_map()` invisibly returns two new fields:
  `zlim` (the color scale used) and `palette` (the palette
  name). Backwards compatible only in the sense that no caller
  was relying on the prior list having exactly four names; tests
  that asserted on `names(...)` were updated.

# rcisignal 0.1.5

## Behavioral change

* `plot_mask_overlay()` now always renders a single panel
  (the base image with the mask overlaid). The
  `side_by_side =` argument is removed (the base-alongside-base
  layout duplicated the base image and added nothing beyond
  what the single overlay panel already shows). The plot region
  is now forced square (`pty = "s"`) to avoid the empty
  canvas that appeared around the image on wide / tall
  devices.

## New features

* New `shift_mask()` exported helper. Slides a logical face-region
  mask by `down` and `right` pixels and returns the shifted mask
  in the same input shape (matrix or column-major vector). This
  was previously a copy-paste recipe in vignette §4.5; it is now
  a real API surface with help page and tests. Useful for
  fine-tuning the *elliptical* `make_face_mask()` regions on a
  non-default base; for the *rectangle* regions, prefer the
  `region_bounds` argument.
* New `region_bounds_from_pixels()` exported helper. Converts
  pixel-coordinate rectangle bounds
  (`c(row_min, row_max, col_min, col_max)` in 1-indexed image
  pixels) to the 0-1 image fractions that
  `make_face_mask(region_bounds = ...)` expects. Bridges visual
  inspection of a base image (where users think in pixels) and
  the parametric API (which expects fractions).
* `plot_mask_overlay()` and `plot_face_mask()` gain a `region =`
  shortcut: `plot_mask_overlay(base, region = "left_eye")`
  builds the mask internally via `make_face_mask()` instead of
  forcing the caller to construct it separately. A
  `region_bounds =` pass-through arg is available for tuning
  the rectangle regions in the same call.

## Tests

* New regression test `test-simulate-diagnostic-chain.R`
  exercises `simulate_2ifc_data()` -> `diagnose_infoval()` /
  `compute_infoval_summary()` / `check_rt_infoval_consistency()`
  / `check_response_inversion()` / `run_diagnostics()` end-to-end
  on the bundled simulator. Guards against recurrence of the
  v0.1.1 rdata-env bugs (CLAUDE.md §11.21: `base_label` and
  missing `noise_type` in the saved rdata file). Skip-gated on
  the `rcicr` Suggests dep and on CRAN-style runs.

## Documentation

* `plot_ci_overlay()` second example switched from
  `signal_region = "eyes"` to `signal_region = "left_eye"` to
  surface the new rectangle eye regions on a frequently-read
  help page (the right-eye rectangle and the rest of the face
  stay untouched, illustrating regional independence).
* Vignette §4.5 mentions of `plot_mask_overlay()` reworded for
  the single-panel design ("overlay on the base image" rather
  than "base alongside base + mask").
* Vignette §4.5 mask-tuning section rewritten to use the new
  exported `shift_mask()` directly instead of redefining it
  inline (the previous copy-paste recipe is gone). The
  `region_bounds` tuning subsection now demos
  `region_bounds_from_pixels()` for the common case where the
  rectangle's edges are eyeballed in pixel coordinates against
  a zoomed-in base image, and points at the
  `plot_mask_overlay(region = ..., region_bounds = ...)`
  shortcut.

## Tooling

* `.gitignore` now excludes `tests/testthat/data`. `rcicr`'s
  `computeInfoVal2IFC()` uses bare `write(...)` for status
  strings, which writes to a file literally named `data` (R's
  default destination for single-argument `write`). The new
  simulate-diagnostic-chain regression test triggers this path
  and would otherwise leave a dirty working tree after every
  test run.

# rcisignal 0.1.4

## Behavioral change

* `make_face_mask(region = "eyes")` now returns a single wide
  axis-aligned **rectangle** covering both eyes ear-to-ear and
  from the eyebrows down to just below the eye line, instead of
  the two small ellipses used in v0.1.0 - v0.1.3. The new
  rectangle is independent of the full-face oval geometry
  (`centre`, `half_width`, `half_height`); rectangle bounds are
  set directly via the new `region_bounds` argument. This makes
  per-eye-line tuning straightforward and removes the previous
  "both eyes move together" coupling that made manual alignment
  to a specific base image awkward. Reverse-correlation results
  computed against the "eyes" region in earlier versions will
  shift numerically because the masked pixel set is broader; the
  Oliveira-2019 worked-example tables in vignette §12 were
  refreshed with the new geometry.

## New features

* New `region = "left_eye"` and `region = "right_eye"` values for
  `make_face_mask()`. Each returns an axis-aligned rectangle
  around the viewer's left or right eye, independent of the
  other eye and of the full-face oval. Bounds default to a
  heuristic centered-face geometry and can be overridden via
  `region_bounds`.
* New `region_bounds = NULL` argument on `make_face_mask()`. Pass
  a length-4 numeric vector `c(row_min, row_max, col_min, col_max)`
  in 0-1 image fractions to override the default bounds for any
  of the three rectangle regions (`"eyes"`, `"left_eye"`,
  `"right_eye"`). Errors if supplied for an elliptical region.

## Documentation

* Vignette §4.5 (Face-region masks) rewritten to introduce the
  eight built-in regions (five ellipses + three rectangles) and
  the two tuning routes: `region_bounds` for rectangles vs
  `centre` / `shift_mask()` for ellipses. The manual-positioning
  example that used to slide both eye ovals via `shift_mask()`
  has been replaced by two examples: a `shift_mask()` demo on
  the mouth ellipse (where the slide workaround is still the
  right tool) and a `region_bounds` demo on the `left_eye`
  rectangle (showing the new precision-tuning knob).
* Vignette face-mask figures regenerated to include `left_eye`
  and `right_eye` panels and the new `mouth_demo_*` /
  `left_eye_demo_*` tuning panels. The old
  `artificial_eyes_demo_*` files were removed.
* README `infoval()` per-region table refreshed: the "Eyes"
  column now reflects the wider rectangle geometry (computed on
  the same Oliveira-2019 Study 1 data, 1000 reference draws,
  matched trial counts). The other columns are unchanged.
* Citation-accuracy sweep: removed several incorrect
  attributions of `make_face_mask()` region geometries to
  Schmitz, Rougier, & Yzerbyt (2024) across the vignette,
  README pointer, and `make_face_mask()` roxygen. Schmitz et
  al. (2024) used an oval mask for `infoval()` calculation but
  did not specify the oval's parameters and did not define any
  sub-region (eyes / nose / mouth / upper-face / lower-face)
  geometry. The package's region defaults are now correctly
  described as this package's own heuristics; the prior
  practice of applying an oval before pixel-wise CI metrics
  is credited to Oliveira et al. (2019), Ratner et al. (2014),
  and Schmitz et al. (2024) jointly. Ratner et al. (2014)
  added to the vignette bibliography. Schmitz et al. (2024) is
  still cited where it applies (Brief-RC structure, `genMask()`
  formula, Brief-RC 12 / 20 variants).

# rcisignal 0.1.3

## New features

* `plot_mask_overlay()` lets you sanity-check that a face mask aligns
  with a specific base image before running CI / reliability /
  discriminability analyses. Accepts a base image (numeric matrix or
  PNG/JPEG path) plus a mask in any of the package's supported forms
  (logical/numeric vector, matrix, or PNG/JPEG path) and draws either
  a side-by-side base + overlay pair or just the overlay panel.
* `plot()` is now defined for `rcisignal_rel_pairwise_report`
  objects returned by `run_discriminability_pairwise()`. One call
  lays out one cluster t-map per pair in a square-ish grid (same
  blue/red sign convention and FWE contours as the per-pair
  `plot.rcisignal_rel_cluster_test()`). A warning fires above
  `max_pairs = 12`; pass `max_pairs = Inf` to silence.
* `plot()` on a `rcisignal_rel_cluster_test` result gains a
  `color_bar = TRUE` argument (default preserves prior behavior).
  Set `FALSE` to suppress the per-panel color bar when packing many
  panels into a small grid; the new pairwise plot method uses this
  internally.

## Documentation

* Vignette §8.4 and §8.5 now show `plot(rep)` as the one-call view
  for the two-condition and pairwise discriminability orchestrators,
  with per-panel `plot(rep$results$...)` shown as the
  custom-title alternative. The `run_discriminability_pairwise()`
  help page is updated in lockstep.

# rcisignal 0.1.2

## New features

* `simulate_2ifc_data()` and `simulate_briefrc_data()` gain an
  `rdata_dir` argument and now return a self-contained `$stimuli`
  list. The simulation object survives `saveRDS()`/`readRDS()`
  and knitr `cache = TRUE` across R sessions: pass an explicit
  `rdata_dir` to keep the stimuli `.Rdata` at a stable path, or
  hand `$stimuli` to downstream consumers in place of
  `$rdata_path`.
* `ci_from_responses_2ifc()`, `diagnose_infoval()`,
  `compute_infoval_summary()`, `check_response_inversion()`,
  `check_rt_infoval_consistency()`, and `run_diagnostics()` gain a
  `stimuli =` argument as an in-memory alternative to
  `rdata_path =` / `rdata =`. When both are supplied `stimuli`
  takes precedence and a warning is emitted.
* `ci_from_responses_briefrc()` argument `base_image_path` is
  renamed `base_image` and now accepts either a numeric matrix
  in `[0, 1]` (e.g. `sim$base_face`) or a file path. The argument
  is now optional when `scaling = "none"` (the default), since
  the base face only feeds the visualization-only `$rendered_ci`
  field. The old `base_image_path` keeps working for one release
  with a deprecation warning.
* `simulate_2ifc_data()` and `simulate_briefrc_data()` now also
  write the base face as a standalone PNG alongside the stimuli
  `.Rdata` (`rcisignal_sim_2ifc_base_face.png` /
  `rcisignal_sim_briefrc_base_face.png` under `rdata_dir`) and
  expose its path on `$base_image_path`.

## Fixes

* `simulate_briefrc_data()` now also writes an rcicr-format
  stimuli `.Rdata` (stable filename
  `rcisignal_sim_briefrc_stimuli.Rdata` under `rdata_dir`) for
  symmetry with the 2IFC path. The Brief-RC consumers still read
  `$noise_matrix` directly so the file is informational rather
  than required.

## Behavioral change

* `plot_ci_overlay()` now uses positive = blue, negative = red
  (matching `plot_agreement_map()` and the cluster-test plots).
  Previously it rendered positive = red, negative = blue, which
  flipped the color reading between the package's three diverging
  plots. Saved overlay PNGs from prior versions will look mirrored
  if regenerated. The four worked-example pairwise figures in
  `vignettes/figures/oliveira_2019/` and the README captions were
  refreshed to match.

## Documentation

* Added `Reading the plot:` sections to `plot_ci_overlay()`,
  `plot_agreement_map()`, `agreement_map_test()`,
  `rel_cluster_test()`, `rel_dissimilarity()`,
  `run_discriminability()`, and `run_discriminability_pairwise()`
  describing color semantics, contour meaning, and
  significance/observed-statistic distinctions.
* Standardized the `mask` `@param` description across
  `agreement_map_test()`, `agreement_map()`, `infoval()`,
  `pixel_t_test()`, `plot_ci_overlay()`, `rel_cluster_test()`,
  `rel_dissimilarity()`, `rel_loo()`, `run_reliability()`,
  `run_discriminability()`, and `run_discriminability_pairwise()`
  to point at `make_face_mask()` / `read_face_mask()` and call
  out the column-major convention.
* Added chained `\dontrun{}` plot examples to `plot_ci_overlay()`
  (with `test = `), `agreement_map_test()` (canonical pairing
  with `plot_ci_overlay()`), `rel_cluster_test()` (cluster t-map
  plot), `rel_dissimilarity()` (bootstrap-distribution plot),
  `run_discriminability()` (whole-report and per-panel plots),
  and `run_discriminability_pairwise()` (per-pair plots and the
  pairwise dissimilarity grid).

# rcisignal 0.1.1

## New features

* `ci_from_responses_briefrc()` now accepts `method = "briefrc20"`
  in addition to `"briefrc12"`. Both Brief-RC variants are
  validated in Schmitz, Rougier, & Yzerbyt (2024). The CI
  computation is identical for both (the `genMask()` formula is
  symmetric in the per-trial split); the `method` argument is
  recorded as `$method` in the result list as provenance metadata.

## Documentation

* Vignette §13.1 rewritten to document Brief-RC 12 and Brief-RC 20
  as the two validated variants. Removed the inaccurate reference
  to a Brief-RC 6 variant (Schmitz et al. mention 4 / 6 / 8 / 10
  only as future research directions, not as published variants).
* Vignette §15 paragraph 5 (group-mean infoVal interpretation)
  rewritten to remove an unsourced "5-10x per-producer median"
  claim and to flag that the sqrt(N)-style inflation of group-mean
  z is conditional on producer-level signal alignment, not
  structural. Brinkman et al. (2019) numerical claims and
  recommendations are retained with explicit page citations.
* Vignette §15 paragraph 1 reframed to drop unsourced
  "10-30% / 70-90%" pixel-fraction figures while keeping the
  Frobenius-norm dilution mechanism.
* README validation section trimmed to a brief pointer; full
  validated-vs-unvalidated breakdown now lives in vignette §1.2.
* New educational background paragraphs added to §7 (split-half /
  Spearman-Brown, ICC), §8 (multiple-comparisons problem, cluster
  permutation, k-connectivity, TFCE, Frobenius/Euclidean
  dissimilarity), and §9 (per-producer and group-mean infoVal).

# rcisignal 0.1.0

First release of the package. Provides a consolidated toolkit for
quality assessment of reverse-correlation data and classification
images: input-side diagnostics (response coding, RT, alignment,
balance) and output-side reliability (pixel-level reliability,
cluster inference, per-producer infoVal).
