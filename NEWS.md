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
