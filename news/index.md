# Changelog

## rcisignal 0.1.4

### Behavioural change

- `make_face_mask(region = "eyes")` now returns a single wide
  axis-aligned **rectangle** covering both eyes ear-to-ear and from the
  eyebrows down to just below the eye line, instead of the two small
  ellipses used in v0.1.0 - v0.1.3. The new rectangle is independent of
  the full-face oval geometry (`centre`, `half_width`, `half_height`);
  rectangle bounds are set directly via the new `region_bounds`
  argument. This makes per-eye-line tuning straightforward and removes
  the previous “both eyes move together” coupling that made manual
  alignment to a specific base image awkward. Reverse-correlation
  results computed against the “eyes” region in earlier versions will
  shift numerically because the masked pixel set is broader; the
  Oliveira-2019 worked-example tables in vignette §12 were refreshed
  with the new geometry.

### New features

- New `region = "left_eye"` and `region = "right_eye"` values for
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md).
  Each returns an axis-aligned rectangle around the viewer’s left or
  right eye, independent of the other eye and of the full-face oval.
  Bounds default to a heuristic centred-face geometry and can be
  overridden via `region_bounds`.
- New `region_bounds = NULL` argument on
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md).
  Pass a length-4 numeric vector `c(row_min, row_max, col_min, col_max)`
  in 0-1 image fractions to override the default bounds for any of the
  three rectangle regions (`"eyes"`, `"left_eye"`, `"right_eye"`).
  Errors if supplied for an elliptical region.

### Documentation

- Vignette §4.5 (Face-region masks) rewritten to introduce the eight
  built-in regions (five ellipses + three rectangles) and the two tuning
  routes: `region_bounds` for rectangles vs `centre` / `shift_mask()`
  for ellipses. The manual-positioning example that used to slide both
  eye ovals via `shift_mask()` has been replaced by two examples: a
  `shift_mask()` demo on the mouth ellipse (where the slide workaround
  is still the right tool) and a `region_bounds` demo on the `left_eye`
  rectangle (showing the new precision-tuning knob).
- Vignette face-mask figures regenerated to include `left_eye` and
  `right_eye` panels and the new `mouth_demo_*` / `left_eye_demo_*`
  tuning panels. The old `artificial_eyes_demo_*` files were removed.
- README
  [`infoval()`](https://olivethree.github.io/rcisignal/reference/infoval.md)
  per-region table refreshed: the “Eyes” column now reflects the wider
  rectangle geometry (computed on the same Oliveira-2019 Study 1 data,
  1000 reference draws, matched trial counts). The other columns are
  unchanged.

## rcisignal 0.1.3

### New features

- [`plot_mask_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_mask_overlay.md)
  lets you sanity-check that a face mask aligns with a specific base
  image before running CI / reliability / discriminability analyses.
  Accepts a base image (numeric matrix or PNG/JPEG path) plus a mask in
  any of the package’s supported forms (logical/numeric vector, matrix,
  or PNG/JPEG path) and draws either a side-by-side base + overlay pair
  or just the overlay panel.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) is now
  defined for `rcisignal_rel_pairwise_report` objects returned by
  [`run_discriminability_pairwise()`](https://olivethree.github.io/rcisignal/reference/run_discriminability_pairwise.md).
  One call lays out one cluster t-map per pair in a square-ish grid
  (same blue/red sign convention and FWE contours as the per-pair
  `plot.rcisignal_rel_cluster_test()`). A warning fires above
  `max_pairs = 12`; pass `max_pairs = Inf` to silence.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  `rcisignal_rel_cluster_test` result gains a `colour_bar = TRUE`
  argument (default preserves prior behaviour). Set `FALSE` to suppress
  the per-panel colour bar when packing many panels into a small grid;
  the new pairwise plot method uses this internally.

### Documentation

- Vignette §8.4 and §8.5 now show `plot(rep)` as the one-call view for
  the two-condition and pairwise discriminability orchestrators, with
  per-panel `plot(rep$results$...)` shown as the custom-title
  alternative. The
  [`run_discriminability_pairwise()`](https://olivethree.github.io/rcisignal/reference/run_discriminability_pairwise.md)
  help page is updated in lockstep.

## rcisignal 0.1.2

### New features

- [`simulate_2ifc_data()`](https://olivethree.github.io/rcisignal/reference/simulate_2ifc_data.md)
  and
  [`simulate_briefrc_data()`](https://olivethree.github.io/rcisignal/reference/simulate_briefrc_data.md)
  gain an `rdata_dir` argument and now return a self-contained
  `$stimuli` list. The simulation object survives
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)/[`readRDS()`](https://rdrr.io/r/base/readRDS.html)
  and knitr `cache = TRUE` across R sessions: pass an explicit
  `rdata_dir` to keep the stimuli `.Rdata` at a stable path, or hand
  `$stimuli` to downstream consumers in place of `$rdata_path`.
- [`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md),
  [`diagnose_infoval()`](https://olivethree.github.io/rcisignal/reference/diagnose_infoval.md),
  [`compute_infoval_summary()`](https://olivethree.github.io/rcisignal/reference/compute_infoval_summary.md),
  [`check_response_inversion()`](https://olivethree.github.io/rcisignal/reference/check_response_inversion.md),
  [`check_rt_infoval_consistency()`](https://olivethree.github.io/rcisignal/reference/check_rt_infoval_consistency.md),
  and
  [`run_diagnostics()`](https://olivethree.github.io/rcisignal/reference/run_diagnostics.md)
  gain a `stimuli =` argument as an in-memory alternative to
  `rdata_path =` / `rdata =`. When both are supplied `stimuli` takes
  precedence and a warning is emitted.
- [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md)
  argument `base_image_path` is renamed `base_image` and now accepts
  either a numeric matrix in `[0, 1]` (e.g. `sim$base_face`) or a file
  path. The argument is now optional when `scaling = "none"` (the
  default), since the base face only feeds the visualisation-only
  `$rendered_ci` field. The old `base_image_path` keeps working for one
  release with a deprecation warning.
- [`simulate_2ifc_data()`](https://olivethree.github.io/rcisignal/reference/simulate_2ifc_data.md)
  and
  [`simulate_briefrc_data()`](https://olivethree.github.io/rcisignal/reference/simulate_briefrc_data.md)
  now also write the base face as a standalone PNG alongside the stimuli
  `.Rdata` (`rcisignal_sim_2ifc_base_face.png` /
  `rcisignal_sim_briefrc_base_face.png` under `rdata_dir`) and expose
  its path on `$base_image_path`.

### Fixes

- [`simulate_briefrc_data()`](https://olivethree.github.io/rcisignal/reference/simulate_briefrc_data.md)
  now also writes an rcicr-format stimuli `.Rdata` (stable filename
  `rcisignal_sim_briefrc_stimuli.Rdata` under `rdata_dir`) for symmetry
  with the 2IFC path. The Brief-RC consumers still read `$noise_matrix`
  directly so the file is informational rather than required.

### Behavioural change

- [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)
  now uses positive = blue, negative = red (matching
  [`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md)
  and the cluster-test plots). Previously it rendered positive = red,
  negative = blue, which flipped the colour reading between the
  package’s three diverging plots. Saved overlay PNGs from prior
  versions will look mirrored if regenerated. The four worked-example
  pairwise figures in `vignettes/figures/oliveira_2019/` and the README
  captions were refreshed to match.

### Documentation

- Added `Reading the plot:` sections to
  [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md),
  [`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md),
  [`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md),
  [`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md),
  [`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md),
  [`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md),
  and
  [`run_discriminability_pairwise()`](https://olivethree.github.io/rcisignal/reference/run_discriminability_pairwise.md)
  describing colour semantics, contour meaning, and
  significance/observed-statistic distinctions.
- Standardised the `mask` `@param` description across
  [`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md),
  `agreement_map()`,
  [`infoval()`](https://olivethree.github.io/rcisignal/reference/infoval.md),
  [`pixel_t_test()`](https://olivethree.github.io/rcisignal/reference/pixel_t_test.md),
  [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md),
  [`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md),
  [`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md),
  [`rel_loo()`](https://olivethree.github.io/rcisignal/reference/rel_loo.md),
  [`run_reliability()`](https://olivethree.github.io/rcisignal/reference/run_reliability.md),
  [`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md),
  and
  [`run_discriminability_pairwise()`](https://olivethree.github.io/rcisignal/reference/run_discriminability_pairwise.md)
  to point at
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  /
  [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
  and call out the column-major convention.
- Added chained `\dontrun{}` plot examples to
  [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)
  (with `test =`),
  [`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md)
  (canonical pairing with
  [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)),
  [`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)
  (cluster t-map plot),
  [`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
  (bootstrap-distribution plot),
  [`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md)
  (whole-report and per-panel plots), and
  [`run_discriminability_pairwise()`](https://olivethree.github.io/rcisignal/reference/run_discriminability_pairwise.md)
  (per-pair plots and the pairwise dissimilarity grid).

## rcisignal 0.1.1

### New features

- [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md)
  now accepts `method = "briefrc20"` in addition to `"briefrc12"`. Both
  Brief-RC variants are validated in Schmitz, Rougier, & Yzerbyt (2024).
  The CI computation is identical for both (the `genMask()` formula is
  symmetric in the per-trial split); the `method` argument is recorded
  as `$method` in the result list as provenance metadata.

### Documentation

- Vignette §13.1 rewritten to document Brief-RC 12 and Brief-RC 20 as
  the two validated variants. Removed the inaccurate reference to a
  Brief-RC 6 variant (Schmitz et al. mention 4 / 6 / 8 / 10 only as
  future research directions, not as published variants).
- Vignette §15 paragraph 5 (group-mean infoVal interpretation) rewritten
  to remove an unsourced “5-10x per-producer median” claim and to flag
  that the sqrt(N)-style inflation of group-mean z is conditional on
  producer-level signal alignment, not structural. Brinkman et
  al. (2019) numerical claims and recommendations are retained with
  explicit page citations.
- Vignette §15 paragraph 1 reframed to drop unsourced “10-30% / 70-90%”
  pixel-fraction figures while keeping the Frobenius-norm dilution
  mechanism.
- README validation section trimmed to a brief pointer; full
  validated-vs-unvalidated breakdown now lives in vignette §1.2.
- New educational background paragraphs added to §7 (split-half /
  Spearman-Brown, ICC), §8 (multiple-comparisons problem, cluster
  permutation, k-connectivity, TFCE, Frobenius/Euclidean dissimilarity),
  and §9 (per-producer and group-mean infoVal).

## rcisignal 0.1.0

First release of the package. Provides a consolidated toolkit for
quality assessment of reverse-correlation data and classification
images: input-side diagnostics (response coding, RT, alignment, balance)
and output-side reliability (pixel-level reliability, cluster inference,
per-producer infoVal).
