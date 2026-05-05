# Simulate 2IFC reverse-correlation data

Generates a synthetic two-image forced-choice (2IFC) dataset that is
shape-compatible with every `check_*()`,
[`run_diagnostics()`](https://olivethree.github.io/rcisignal/reference/run_diagnostics.md),
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md),
and reliability/discriminability function in the package. Useful as a
quickstart sandbox and as a building block for simulation studies
(power, calibration of reliability and discriminability metrics,
sensitivity to contamination).

The function generates the noise pool on the fly via
[`rcicr::generateNoisePattern()`](https://rdrr.io/pkg/rcicr/man/generateNoisePattern.html)
and
[`rcicr::generateNoiseImage()`](https://rdrr.io/pkg/rcicr/man/generateNoiseImage.html),
so it requires the `rcicr` package to be installed. With default
arguments (50 participants per condition, 500 trials, 256-pixel images)
the noise generation step takes one to a few minutes; a progress bar is
shown.

## Usage

``` r
simulate_2ifc_data(
  n_per_condition = 50L,
  conditions = c("target", "control"),
  n_trials = 500L,
  img_size = 256L,
  base_image = NULL,
  signal_strength = "weak",
  signal_region = "eyes",
  rt_contamination_fast = 0.02,
  rt_contamination_slow = 0.02,
  noise_type = "sinusoid",
  nscales = 5L,
  sigma = 25,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- n_per_condition:

  Integer. Participants per condition. Default `50`.

- conditions:

  Character vector. Default `c("target", "control")`.

- n_trials:

  Integer. Trials per participant; equals the noise pool size (each pool
  item shown once per participant). Default `500`.

- img_size:

  Integer. Side length of square images, in pixels. Default `256`
  (matches the bundled base face). Setting this higher (e.g. `512`)
  requires you to also pass a matching `base_image`.

- base_image:

  Path to a square PNG, or a numeric matrix in `[0, 1]` of dimension
  `img_size x img_size`. Default `NULL` uses the bundled
  `inst/extdata/sim_base_face.png` (a 256x256 grayscale artificial
  face).

- signal_strength:

  One of `"none"`, `"weak"`, `"strong"`, or a numeric coefficient (the
  `beta` in the logistic above). Default `"weak"`.

- signal_region:

  Region passed to
  [`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md).
  Default `"eyes"`.

- rt_contamination_fast, rt_contamination_slow:

  Numeric in `[0, 1]`. Fraction of trials replaced by uniform-fast
  (50-200 ms) and uniform-slow (5000-20000 ms) contaminants. Default
  `0.02` each.

- noise_type, nscales, sigma:

  Forwarded to
  [`rcicr::generateNoisePattern()`](https://rdrr.io/pkg/rcicr/man/generateNoisePattern.html).
  Defaults match rcicr's defaults (`"sinusoid"`, `5`, `25`).

- seed:

  Integer or `NULL`. When `NULL`, a random seed is drawn and stored on
  the result so the run is reproducible.

- progress:

  Logical. Show a `cli` progress bar during noise generation. Default
  `TRUE`.

## Value

An object of class `"rcisignal_sim"` with elements:

- `data` — a
  [data.table::data.table](https://rdrr.io/pkg/data.table/man/data.table.html)
  with one row per trial and columns `participant_id`, `condition`,
  `trial`, `stimulus`, `response`, `rt`. Compatible with
  [`run_diagnostics()`](https://olivethree.github.io/rcisignal/reference/run_diagnostics.md)
  and
  [`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md).

- `noise_matrix` — `(img_size^2) x n_trials` numeric matrix.

- `base_face` — `img_size x img_size` numeric matrix.

- `params` — `n_trials x ncoef` matrix of per-trial sinusoid
  coefficients (the rcicr `stimuli_params`).

- `p` — the rcicr noise basis (`generateNoisePattern()` output); pair
  with `params` to regenerate any noise image.

- `rdata_path` — path to an rcicr-format `.Rdata` file written to a
  session tempdir, suitable for
  [`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md)
  /
  [`compute_infoval_summary()`](https://olivethree.github.io/rcisignal/reference/compute_infoval_summary.md)
  and other downstream functions that take an `rdata` argument.

- `signal` — pixel-level signal vector used to plant the response bias.

- `meta` — list of method, `n_per_condition`, `conditions`, `n_trials`,
  `img_size`, `signal_strength`, `signal_region`, `seed`,
  `elapsed_secs`.

## Signal model

Each trial `t` shows two stimuli, image_a = base +
noise[t](https://rdrr.io/r/base/t.html) and image_b = base -
noise[t](https://rdrr.io/r/base/t.html). The participant chooses one.
With `signal_strength = "none"` choices are uniform random
(`P(+1) = 0.5`); with `"weak"` / `"strong"` (or a custom numeric
`beta`), the log-odds of choosing image_a are
`beta * (noise[, t] %*% s) / scale`, where `s` is a binary mask from
[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md)
over the requested `signal_region` and `scale = sqrt(sum(s))` keeps the
per-pixel signal magnitude comparable across regions of different size.

## RT model

Shifted lognormal: `rt = round(exp(rnorm(n, log(800), 0.5)) + 150)`, in
ms. A small fraction of fast (\<200 ms) and slow (\>5000 ms)
contaminants are mixed in (default 2% each) so the diagnostic functions
([`check_rt()`](https://olivethree.github.io/rcisignal/reference/check_rt.md))
have something to flag.

## See also

[`simulate_briefrc_data()`](https://olivethree.github.io/rcisignal/reference/simulate_briefrc_data.md),
[`run_diagnostics()`](https://olivethree.github.io/rcisignal/reference/run_diagnostics.md),
[`ci_from_responses_2ifc()`](https://olivethree.github.io/rcisignal/reference/ci_from_responses_2ifc.md),
[`run_reliability()`](https://olivethree.github.io/rcisignal/reference/run_reliability.md),
[`run_discriminability_pairwise()`](https://olivethree.github.io/rcisignal/reference/run_discriminability_pairwise.md),
[`make_face_mask()`](https://olivethree.github.io/rcisignal/reference/make_face_mask.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# `sim$data` is a plain data frame (columns: participant_id, stimulus,
# response, condition, rt) — same shape ci_from_responses_2ifc() and
# the check_*() functions expect from your own CSV.
sim <- simulate_2ifc_data(n_per_condition = 10, n_trials = 60, seed = 1)
run_diagnostics(sim$data, method = "2ifc", col_rt = "rt")
cis <- ci_from_responses_2ifc(sim$data, rdata_path = sim$rdata_path)
run_reliability(cis$signal_matrix, n_permutations = 200L, seed = 1)
} # }
```
