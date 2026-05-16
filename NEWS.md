# rcisignal 0.1.6

## New features

* `plot_agreement_map()` gains a unipolar `palette = "fire"` option
  that displays `|t|` on a single-hue (pale yellow to deep red)
  ramp. Use when the question is "where do producers have a
  consistent opinion" and direction is not needed; the `"fire"`
  view discards sign by design. The default `palette = "diverging"`
  is unchanged.

## Breaking change

* `plot_agreement_map(palette = "viridis")` is removed. The
  previous viridis branch plotted signed `t` against a palette
  with no neutral midpoint, so zero rendered at a perceptually
  arbitrary colour. The new `"fire"` option supersedes that use
  case correctly. The argument is now
  `palette = c("diverging", "fire")`.

## Documentation

* `plot_agreement_map()` help page rewritten to be explicit that
  on the default diverging palette both deep red and deep blue
  indicate **strong** agreement (sign in hue, magnitude in
  saturation); "no agreement" is the neutral colour (white), not
  red. The new `"fire"` palette is documented as a sign-discarding
  magnitude-only view with explicit pointers to recover direction
  via the diverging palette or `plot_ci_overlay()`.
* Vignette §10.2 gains a one-paragraph aside introducing the
  fire option and reiterating the two-channel reading of the
  diverging palette.

## Internal

* `plot_agreement_map()` invisibly returns two new fields:
  `zlim` (the colour scale used) and `palette` (the palette
  name). Backwards compatible only in the sense that no caller
  was relying on the prior list having exactly four names; tests
  that asserted on `names(...)` were updated.

# rcisignal 0.1.5

## Behavioural change

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

## Behavioural change

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
  heuristic centred-face geometry and can be overridden via
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
  `colour_bar = TRUE` argument (default preserves prior behaviour).
  Set `FALSE` to suppress the per-panel colour bar when packing many
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
  the base face only feeds the visualisation-only `$rendered_ci`
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

## Behavioural change

* `plot_ci_overlay()` now uses positive = blue, negative = red
  (matching `plot_agreement_map()` and the cluster-test plots).
  Previously it rendered positive = red, negative = blue, which
  flipped the colour reading between the package's three diverging
  plots. Saved overlay PNGs from prior versions will look mirrored
  if regenerated. The four worked-example pairwise figures in
  `vignettes/figures/oliveira_2019/` and the README captions were
  refreshed to match.

## Documentation

* Added `Reading the plot:` sections to `plot_ci_overlay()`,
  `plot_agreement_map()`, `agreement_map_test()`,
  `rel_cluster_test()`, `rel_dissimilarity()`,
  `run_discriminability()`, and `run_discriminability_pairwise()`
  describing colour semantics, contour meaning, and
  significance/observed-statistic distinctions.
* Standardised the `mask` `@param` description across
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
