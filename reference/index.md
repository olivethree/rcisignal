# Package index

## Build per-producer CIs

Compute one classification image (CI) per producer from trial-level
responses. Pass `group_by =` to also build group-averaged CIs in one
call.

- [`ci_from_responses_briefrc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_briefrc.md)
  : Compute individual Brief-RC CIs from trial-level responses
- [`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md)
  : Compute individual 2IFC CIs from trial-level responses
- [`group_ci()`](https://olivethree.github.io/rcisignal/reference/group_ci.md)
  : Collapse a per-producer signal matrix into per-group means

## Reliability and informational value

Within-condition reliability of the per-producer signal and per-producer
infoVal.

- [`infoval()`](https://olivethree.github.io/rcisignal/reference/infoval.md)
  : Per-producer informational value with trial-count-matched reference
- [`rel_split_half()`](https://olivethree.github.io/rcisignal/reference/rel_split_half.md)
  : Permuted split-half reliability with Spearman-Brown correction
- [`rel_split_half_null()`](https://olivethree.github.io/rcisignal/reference/rel_split_half_null.md)
  : Compute only the empirical-null distribution for split-half
- [`rel_icc()`](https://olivethree.github.io/rcisignal/reference/rel_icc.md)
  : Intraclass correlation coefficients via direct mean squares
- [`rel_loo()`](https://olivethree.github.io/rcisignal/reference/rel_loo.md)
  : Leave-one-out influence screening
- [`rel_loo_z()`](https://olivethree.github.io/rcisignal/reference/rel_loo_z.md)
  : Z-scored leave-one-out influence (accessor)
- [`run_reliability()`](https://olivethree.github.io/rcisignal/reference/run_reliability.md)
  : Run every within-condition reliability metric

## Between-condition discriminability

Between-condition spatial difference: pixel-wise t, cluster-based
permutation, Frobenius-distance bootstrap.

- [`pixel_t_test()`](https://olivethree.github.io/rcisignal/reference/pixel_t_test.md)
  : Vectorised pixel-wise t-test (independent or paired)
- [`rel_cluster_test()`](https://olivethree.github.io/rcisignal/reference/rel_cluster_test.md)
  : Cluster-based permutation test with family-wise error control
- [`rel_dissimilarity()`](https://olivethree.github.io/rcisignal/reference/rel_dissimilarity.md)
  : Between-condition dissimilarity with bootstrap confidence intervals
- [`run_discriminability()`](https://olivethree.github.io/rcisignal/reference/run_discriminability.md)
  : Run every between-condition discriminability metric
- [`run_discriminability_pairwise()`](https://olivethree.github.io/rcisignal/reference/run_discriminability_pairwise.md)
  : Run all pairwise between-condition comparisons across K conditions

## Cross-producer agreement

Producer-level agreement maps (descriptive and FWE-controlled).

- [`agreement_map_test()`](https://olivethree.github.io/rcisignal/reference/agreement_map_test.md)
  : Per-pixel inferential agreement map with FWER control
- [`plot_agreement_map()`](https://olivethree.github.io/rcisignal/reference/plot_agreement_map.md)
  : Per-pixel agreement heatmap for a producer signal matrix

## Compare multiple CIs

Cross-CI similarity, distance, and projection views. Accept any named
collection of CIs (per-producer or group-level).

- [`plot_ci_correlogram()`](https://olivethree.github.io/rcisignal/reference/plot_ci_correlogram.md)
  : Pairwise Pearson-r correlogram across multiple group-mean CIs
- [`plot_ci_distance_matrix()`](https://olivethree.github.io/rcisignal/reference/plot_ci_distance_matrix.md)
  : Pairwise Euclidean distance matrix across multiple group-mean CIs
- [`plot_ci_mds()`](https://olivethree.github.io/rcisignal/reference/plot_ci_mds.md)
  : Multidimensional-scaling (MDS) projection of multiple CIs

## Plot and export CIs

- [`plot_ci_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_ci_overlay.md)
  : Render a group CI as a translucent agreement-map overlay on a base
  image
- [`plot_mask_overlay()`](https://olivethree.github.io/rcisignal/reference/plot_mask_overlay.md)
  : Sanity-check a face mask against a base image
- [`plot_face_mask()`](https://olivethree.github.io/rcisignal/reference/plot_face_mask.md)
  : Plot a face mask for visual verification
- [`plot_dissimilarity_grid()`](https://olivethree.github.io/rcisignal/reference/plot_dissimilarity_grid.md)
  : Compare Euclidean dissimilarity intervals across multiple contrasts
  on a single plot
- [`save_ci_images()`](https://olivethree.github.io/rcisignal/reference/save_ci_images.md)
  : Write rendered CIs to PNG or JPEG files

## I/O and simulation

- [`read_responses()`](https://olivethree.github.io/rcisignal/reference/read_responses.md)
  : Read reverse-correlation response data from a CSV file
- [`read_noise_matrix()`](https://olivethree.github.io/rcisignal/reference/read_noise_matrix.md)
  : Read a noise matrix from any supported source, with caching
- [`read_cis()`](https://olivethree.github.io/rcisignal/reference/read_cis.md)
  : Read a directory of CI images into a pixels x participants matrix
- [`read_signal_matrix()`](https://olivethree.github.io/rcisignal/reference/read_signal_matrix.md)
  : Read a directory of CI images and extract the base-subtracted signal
- [`extract_signal()`](https://olivethree.github.io/rcisignal/reference/extract_signal.md)
  : Subtract the base image from each column of a CI matrix
- [`simulate_briefrc_data()`](https://olivethree.github.io/rcisignal/reference/simulate_briefrc_data.md)
  : Simulate Brief-RC reverse-correlation data
- [`simulate_2ifc_data()`](https://olivethree.github.io/rcisignal/reference/simulate_2ifc_data.md)
  : Simulate 2IFC reverse-correlation data

## Data-quality diagnostics

- [`run_diagnostics()`](https://olivethree.github.io/rcisignal/reference/run_diagnostics.md)
  : Run the full battery of diagnostic checks
- [`infoval_report()`](https://olivethree.github.io/rcisignal/reference/infoval_report.md)
  : Per-producer infoVal report (with calibration sanity checks)
- [`check_duplicates()`](https://olivethree.github.io/rcisignal/reference/check_duplicates.md)
  : Check for duplicate rows in response data
- [`check_response_bias()`](https://olivethree.github.io/rcisignal/reference/check_response_bias.md)
  : Check for response bias
- [`check_response_coding()`](https://olivethree.github.io/rcisignal/reference/check_response_coding.md)
  : Check that response values use the expected coding
- [`check_response_inversion()`](https://olivethree.github.io/rcisignal/reference/check_response_inversion.md)
  : Detect response-code inversion
- [`check_rt()`](https://olivethree.github.io/rcisignal/reference/check_rt.md)
  : Response-time distributional checks
- [`check_rt_infoval_consistency()`](https://olivethree.github.io/rcisignal/reference/check_rt_infoval_consistency.md)
  : Cross-validate infoVal against RT quality
- [`check_stimulus_alignment()`](https://olivethree.github.io/rcisignal/reference/check_stimulus_alignment.md)
  : Check stimulus ids align with the stimulus pool
- [`check_trial_counts()`](https://olivethree.github.io/rcisignal/reference/check_trial_counts.md)
  : Check trial counts per participant
- [`check_version_compat()`](https://olivethree.github.io/rcisignal/reference/check_version_compat.md)
  : Check rcicr version compatibility with a 2IFC rdata file
- [`validate_noise_matrix()`](https://olivethree.github.io/rcisignal/reference/validate_noise_matrix.md)
  : Validate a noise matrix

## Helpers

- [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
  : Build an oval face-region mask for a square image

- [`read_face_mask()`](https://olivethree.github.io/rcisignal/reference/read_face_mask.md)
  : Read an image-based face-region mask from a PNG or JPEG

- [`shift_mask()`](https://olivethree.github.io/rcisignal/reference/shift_mask.md)
  : Shift a face-region mask by a number of pixels

- [`region_bounds_from_pixels()`](https://olivethree.github.io/rcisignal/reference/region_bounds_from_pixels.md)
  : Convert pixel-coordinate rectangle bounds to 0-1 image fractions

- [`is_rcisignal_diag_result()`](https://olivethree.github.io/rcisignal/reference/is_rcisignal_diag_result.md)
  :

  Test whether an object is an `rcisignal_diag_result`

- [`rcisignal_diag_result()`](https://olivethree.github.io/rcisignal/reference/rcisignal_diag_result.md)
  :

  Construct an `rcisignal_diag_result` object

- [`rcisignal_diag_report()`](https://olivethree.github.io/rcisignal/reference/rcisignal_diag_report.md)
  :

  Construct an `rcisignal_diag_report` object
