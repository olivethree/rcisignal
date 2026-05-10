#' Simulate 2IFC reverse-correlation data
#'
#' @description
#' Generates a synthetic two-image forced-choice (2IFC) dataset that
#' is shape-compatible with every `check_*()`, `run_diagnostics()`,
#' `ci_from_responses_2ifc()`, and reliability/discriminability
#' function in the package. Useful as a quickstart sandbox and as a
#' building block for simulation studies (power, calibration of
#' reliability and discriminability metrics, sensitivity to
#' contamination).
#'
#' The function generates the noise pool on the fly via
#' [rcicr::generateNoisePattern()] and [rcicr::generateNoiseImage()],
#' so it requires the `rcicr` package to be installed. With default
#' arguments (50 participants per condition, 500 trials, 256-pixel
#' images) the noise generation step takes one to a few minutes; a
#' progress bar is shown.
#'
#' @section Signal model:
#' Each trial `t` shows two stimuli, image_a = base + noise[t] and
#' image_b = base - noise[t]. The participant chooses one. With
#' `signal_strength = "none"` choices are uniform random (`P(+1) =
#' 0.5`); with `"weak"` / `"strong"` (or a custom numeric `beta`),
#' the log-odds of choosing image_a are
#' `beta * (noise[, t] %*% s) / scale`, where `s` is a binary mask
#' from [make_face_mask()] over the requested `signal_region` and
#' `scale = sqrt(sum(s))` keeps the per-pixel signal magnitude
#' comparable across regions of different size.
#'
#' @section RT model:
#' Shifted lognormal: `rt = round(exp(rnorm(n, log(800), 0.5)) +
#' 150)`, in ms. A small fraction of fast (<200 ms) and slow
#' (>5000 ms) contaminants are mixed in (default 2% each) so the
#' diagnostic functions ([check_rt()]) have something to flag.
#'
#' @param n_per_condition Integer. Participants per condition.
#'   Default `50`.
#' @param conditions Character vector. Default
#'   `c("target", "control")`.
#' @param n_trials Integer. Trials per participant; equals the noise
#'   pool size (each pool item shown once per participant). Default
#'   `500`.
#' @param img_size Integer. Side length of square images, in pixels.
#'   Default `256` (matches the bundled base face). Setting this
#'   higher (e.g. `512`) requires you to also pass a matching
#'   `base_image`.
#' @param base_image Path to a square PNG, or a numeric matrix in
#'   `[0, 1]` of dimension `img_size x img_size`. Default `NULL`
#'   uses the bundled `inst/extdata/sim_base_face.png` (a 256x256
#'   grayscale artificial face).
#' @param signal_strength One of `"none"`, `"weak"`, `"strong"`, or
#'   a numeric coefficient (the `beta` in the logistic above).
#'   Default `"weak"`.
#' @param signal_region Region passed to [make_face_mask()]. Default
#'   `"eyes"`.
#' @param rt_contamination_fast,rt_contamination_slow Numeric in
#'   `[0, 1]`. Fraction of trials replaced by uniform-fast (50-200
#'   ms) and uniform-slow (5000-20000 ms) contaminants. Default
#'   `0.02` each.
#' @param noise_type,nscales,sigma Forwarded to
#'   [rcicr::generateNoisePattern()]. Defaults match rcicr's
#'   defaults (`"sinusoid"`, `5`, `25`).
#' @param rdata_dir Optional directory in which to write the
#'   rcicr-format `.Rdata` stimuli file with a stable filename
#'   (`rcisignal_sim_2ifc_stimuli.Rdata`). When `NULL` (default)
#'   the file goes to a session tempdir and the returned
#'   `$rdata_path` becomes invalid after the R session ends. Pass
#'   an explicit directory to persist the simulation across
#'   sessions. See Details.
#' @param seed Integer or `NULL`. When `NULL`, a random seed is
#'   drawn and stored on the result so the run is reproducible.
#' @param progress Logical. Show a `cli` progress bar during noise
#'   generation. Default `TRUE`.
#'
#' @details
#' When `rdata_dir = NULL`, the returned `$rdata_path` points at a
#' session tempdir and becomes invalid after the R session ends.
#' Persist the simulation across sessions (caching with
#' [saveRDS()], knitr `cache=TRUE`, sharing with collaborators)
#' either by passing an explicit `rdata_dir` such as `"simdata/"`,
#' or by handing the returned `$stimuli` list to
#' [ci_from_responses_2ifc()] (and the other infoval-dependent
#' helpers) in place of `rdata_path`. `$stimuli` is a
#' self-contained in-memory representation of the rcicr stimuli
#' env, so the sim object round-trips through
#' `saveRDS()`/`readRDS()` without a file dependency.
#'
#' @return An object of class `"rcisignal_sim"` with elements:
#'   * `data` — a [data.table::data.table] with one row per trial
#'     and columns `participant_id`, `condition`, `trial`,
#'     `stimulus`, `response`, `rt`. Compatible with
#'     [run_diagnostics()] and [ci_from_responses_2ifc()].
#'   * `noise_matrix` — `(img_size^2) x n_trials` numeric matrix.
#'   * `base_face` — `img_size x img_size` numeric matrix.
#'   * `params` — `n_trials x ncoef` matrix of per-trial sinusoid
#'     coefficients (the rcicr `stimuli_params`).
#'   * `p` — the rcicr noise basis (`generateNoisePattern()`
#'     output); pair with `params` to regenerate any noise image.
#'   * `rdata_path` — path to an rcicr-format `.Rdata` file
#'     written either to a session tempdir (when
#'     `rdata_dir = NULL`) or to the user-supplied `rdata_dir`.
#'     Suitable for [ci_from_responses_2ifc()] /
#'     [compute_infoval_summary()] and other downstream functions
#'     that take an `rdata` argument. **Not portable across R
#'     sessions** when `rdata_dir = NULL`.
#'   * `stimuli` — a self-contained list (`base_face`, `params`,
#'     `p`, etc.) that downstream consumers accept via their
#'     `stimuli =` argument as a portable alternative to
#'     `rdata_path`. Round-trips through
#'     [saveRDS()]/[readRDS()].
#'   * `signal` — pixel-level signal vector used to plant the
#'     response bias.
#'   * `meta` — list of method, `n_per_condition`, `conditions`,
#'     `n_trials`, `img_size`, `signal_strength`, `signal_region`,
#'     `seed`, `elapsed_secs`.
#'
#' @seealso [simulate_briefrc_data()], [run_diagnostics()],
#'   [ci_from_responses_2ifc()], [run_reliability()],
#'   [run_discriminability_pairwise()], [make_face_mask()].
#'
#' @examples
#' \dontrun{
#' # `sim$data` is a plain data frame (columns: participant_id, stimulus,
#' # response, condition, rt) — same shape ci_from_responses_2ifc() and
#' # the check_*() functions expect from your own CSV.
#' sim <- simulate_2ifc_data(n_per_condition = 10, n_trials = 60, seed = 1)
#' run_diagnostics(sim$data, method = "2ifc", col_rt = "rt")
#' cis <- ci_from_responses_2ifc(sim$data, rdata_path = sim$rdata_path)
#' run_reliability(cis$signal_matrix, n_permutations = 200L, seed = 1)
#' }
#' @export
simulate_2ifc_data <- function(n_per_condition       = 50L,
                               conditions            = c("target",
                                                         "control"),
                               n_trials              = 500L,
                               img_size              = 256L,
                               base_image            = NULL,
                               signal_strength       = "weak",
                               signal_region         = "eyes",
                               rt_contamination_fast = 0.02,
                               rt_contamination_slow = 0.02,
                               noise_type            = "sinusoid",
                               nscales               = 5L,
                               sigma                 = 25,
                               rdata_dir             = NULL,
                               seed                  = NULL,
                               progress              = TRUE) {
  ensure_rcicr()
  validate_simulate_args(
    n_per_condition, conditions, n_trials, img_size,
    rt_contamination_fast, rt_contamination_slow
  )
  beta <- resolve_signal_strength(signal_strength)

  used_seed <- resolve_seed(seed)
  old_seed  <- save_rng_state()
  on.exit(restore_rng_state(old_seed), add = TRUE)
  set.seed(used_seed)

  t0 <- Sys.time()
  base_face <- load_base_face(base_image, img_size)
  pool <- generate_noise_pool(
    n_trials = n_trials, img_size = img_size,
    noise_type = noise_type, nscales = nscales, sigma = sigma,
    progress = progress
  )

  signal <- make_signal_vector(img_size, signal_region)
  scale  <- sqrt(sum(signal))
  if (scale == 0) scale <- 1

  data_list <- vector("list", length(conditions) * n_per_condition)
  k <- 1L
  for (cond in conditions) {
    for (pid in seq_len(n_per_condition)) {
      logits <- if (beta == 0) rep(0, n_trials) else
        beta * as.vector(crossprod(pool$noise_matrix, signal)) / scale
      probs   <- stats::plogis(logits)
      choices <- stats::rbinom(n_trials, 1L, probs)
      response <- ifelse(choices == 1L, 1L, -1L)
      rts <- simulate_rts(
        n_trials,
        rt_contamination_fast,
        rt_contamination_slow
      )
      data_list[[k]] <- data.table::data.table(
        participant_id = sprintf("%s_%03d", cond, pid),
        condition      = cond,
        trial          = seq_len(n_trials),
        stimulus       = seq_len(n_trials),
        response       = response,
        rt             = rts
      )
      k <- k + 1L
    }
  }
  data <- data.table::rbindlist(data_list)

  stimuli <- build_sim_stimuli(
    base_face  = base_face, p = pool$p, params = pool$params,
    img_size   = img_size, n_trials = n_trials, seed = used_seed,
    noise_type = noise_type, nscales = nscales, sigma = sigma
  )
  rdata_path <- write_sim_rdata(
    stimuli = stimuli, dir = rdata_dir, method = "2ifc"
  )
  inform_rdata_written(rdata_path, persistent = !is.null(rdata_dir))

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  new_rcisignal_sim(
    data         = data,
    noise_matrix = pool$noise_matrix,
    base_face    = base_face,
    params       = pool$params,
    p            = pool$p,
    signal       = signal,
    rdata_path   = rdata_path,
    stimuli      = stimuli,
    meta         = list(
      method            = "2ifc",
      n_per_condition   = as.integer(n_per_condition),
      conditions        = as.character(conditions),
      n_trials          = as.integer(n_trials),
      img_size          = as.integer(img_size),
      signal_strength   = signal_strength,
      signal_region     = signal_region,
      seed              = used_seed,
      elapsed_secs      = elapsed
    )
  )
}

#' Simulate Brief-RC reverse-correlation data
#'
#' @description
#' Generates a synthetic Brief-RC dataset (Schmitz, Rougier, &
#' Yzerbyt, 2024) where each trial shows several original/inverted
#' noise pairs and the participant picks one image. Output is
#' shape-compatible with [run_diagnostics()],
#' [ci_from_responses_briefrc()], and the reliability /
#' discriminability functions.
#'
#' Noise pool is generated once via [rcicr::generateNoisePattern()]
#' and [rcicr::generateNoiseImage()] and then sampled per trial
#' (without replacement within a participant; the same pool is
#' shared across participants).
#'
#' @section Signal model:
#' Per trial, with `images_per_trial = 2k`, the participant sees `k`
#' original/inverted pairs. Each of the `2k` images has utility
#' `beta * (noise %*% s) / scale + Gumbel(0, 1)`, where `noise` is
#' `+noise[, j]` for the original version and `-noise[, j]` for the
#' inverted version of pair `j`. The participant picks the image
#' with the highest utility (multinomial logit / softmax). The
#' recorded `stimulus` is the pool index of the chosen pair;
#' `response` is `+1` if the original version was chosen and `-1`
#' for the inverted version.
#'
#' @inheritSection simulate_2ifc_data RT model
#'
#' @param n_per_condition,conditions,img_size,base_image,signal_strength,signal_region,rt_contamination_fast,rt_contamination_slow,noise_type,nscales,sigma,rdata_dir,seed,progress
#'   See [simulate_2ifc_data()]. For Brief-RC the rcicr-format
#'   `.Rdata` is informational (downstream Brief-RC functions read
#'   `$noise_matrix` directly); it is written for symmetry with the
#'   2IFC path and as a portable on-disk artefact when
#'   `rdata_dir` is non-`NULL` (stable filename
#'   `rcisignal_sim_briefrc_stimuli.Rdata`).
#' @param n_trials Integer or `NULL`. Brief-RC trials per
#'   participant. When `NULL` (default), it is derived from the
#'   pair budget as `noise_pool_size %/% (images_per_trial / 2)`,
#'   so the simulation matches `simulate_2ifc_data()` on number of
#'   unique noise pairs rather than on trial count. With the
#'   defaults (`noise_pool_size = 500`, `images_per_trial = 12`)
#'   this resolves to `n_trials = 83`; with `images_per_trial =
#'   20` it resolves to `50`. Supplying `n_trials` directly
#'   overrides this and (if `noise_pool_size` is also `NULL`)
#'   restores the older "pool grows with trials" behaviour.
#' @param images_per_trial Integer (even). Number of images shown
#'   per trial; half are original and half are inverted versions of
#'   the same noise patterns. Default `12` (= 6 pairs).
#' @param noise_pool_size Integer or `NULL`. Total number of noise
#'   patterns to pre-generate. When `NULL` and `n_trials` is also
#'   `NULL`, defaults to `500` (matched image-pair budget against
#'   the 2IFC default of 500 trials x 1 pair). When `NULL` but
#'   `n_trials` is supplied, defaults to `n_trials *
#'   (images_per_trial / 2)` so within-participant sampling stays
#'   without replacement. Pass a larger value than that to study
#'   sub-sampling.
#'
#' @return An object of class `"rcisignal_sim"`. See
#'   [simulate_2ifc_data()] for the structure. The `meta` list also
#'   carries `images_per_trial` and `noise_pool_size`.
#'
#' @references
#' Schmitz, M., Rougier, M., & Yzerbyt, V. (2024). Introducing the
#' brief reverse correlation: an improved tool to assess visual
#' representations. *European Journal of Social Psychology*.
#' \doi{10.1002/ejsp.3100}
#'
#' @seealso [simulate_2ifc_data()], [run_diagnostics()],
#'   [ci_from_responses_briefrc()].
#'
#' @examples
#' \dontrun{
#' # `sim$data` is a plain data frame (columns: participant_id, stimulus,
#' # response, condition, rt) — same shape ci_from_responses_briefrc() and
#' # the check_*() functions expect from your own CSV.
#' # `sim$noise_matrix` is a numeric matrix `n_pixels x pool_size` — same
#' # shape read_noise_matrix() returns from your Brief-RC OSF txt file.
#' sim <- simulate_briefrc_data(n_per_condition = 10, n_trials = 60, seed = 1)
#' run_diagnostics(sim$data, method = "briefrc",
#'                 noise_matrix = sim$noise_matrix, col_rt = "rt")
#' cis <- ci_from_responses_briefrc(sim$data,
#'                                  noise_matrix = sim$noise_matrix)
#' run_reliability(cis$signal_matrix, n_permutations = 200L, seed = 1)
#' }
#' @export
simulate_briefrc_data <- function(n_per_condition       = 50L,
                                  conditions            = c("target",
                                                            "control"),
                                  n_trials              = NULL,
                                  images_per_trial      = 12L,
                                  noise_pool_size       = NULL,
                                  img_size              = 256L,
                                  base_image            = NULL,
                                  signal_strength       = "weak",
                                  signal_region         = "eyes",
                                  rt_contamination_fast = 0.02,
                                  rt_contamination_slow = 0.02,
                                  noise_type            = "sinusoid",
                                  nscales               = 5L,
                                  sigma                 = 25,
                                  rdata_dir             = NULL,
                                  seed                  = NULL,
                                  progress              = TRUE) {
  ensure_rcicr()
  if (!is.numeric(images_per_trial) || length(images_per_trial) != 1L ||
        images_per_trial < 2L || images_per_trial %% 2L != 0L) {
    cli::cli_abort(
      "{.arg images_per_trial} must be an even integer >= 2."
    )
  }
  pairs_per_trial <- as.integer(images_per_trial) %/% 2L

  if (is.null(n_trials) && is.null(noise_pool_size)) {
    noise_pool_size <- 500L
    n_trials <- noise_pool_size %/% pairs_per_trial
  } else if (is.null(n_trials)) {
    if (!is.numeric(noise_pool_size) || length(noise_pool_size) != 1L) {
      cli::cli_abort(
        "{.arg noise_pool_size} must be a single integer or NULL."
      )
    }
    n_trials <- as.integer(noise_pool_size) %/% pairs_per_trial
  } else if (is.null(noise_pool_size)) {
    noise_pool_size <- as.integer(n_trials) * pairs_per_trial
  }

  validate_simulate_args(
    n_per_condition, conditions, n_trials, img_size,
    rt_contamination_fast, rt_contamination_slow
  )
  if (!is.numeric(noise_pool_size) || length(noise_pool_size) != 1L ||
        noise_pool_size < as.integer(n_trials) * pairs_per_trial) {
    cli::cli_abort(c(
      "{.arg noise_pool_size} must be a single integer \\
       >= {.val {as.integer(n_trials) * pairs_per_trial}} \\
       (= {.arg n_trials} x {.arg images_per_trial} / 2).",
      "i" = "Within-participant sampling is without replacement, \\
             so the pool must hold at least one full session's \\
             worth of pairs."
    ))
  }
  beta <- resolve_signal_strength(signal_strength)

  used_seed <- resolve_seed(seed)
  old_seed  <- save_rng_state()
  on.exit(restore_rng_state(old_seed), add = TRUE)
  set.seed(used_seed)

  t0 <- Sys.time()
  base_face <- load_base_face(base_image, img_size)
  pool <- generate_noise_pool(
    n_trials   = noise_pool_size, img_size = img_size,
    noise_type = noise_type, nscales = nscales, sigma = sigma,
    progress   = progress
  )

  signal <- make_signal_vector(img_size, signal_region)
  scale  <- sqrt(sum(signal))
  if (scale == 0) scale <- 1

  proj <- as.vector(crossprod(pool$noise_matrix, signal)) / scale

  data_list <- vector("list", length(conditions) * n_per_condition)
  k <- 1L
  draw_n <- as.integer(n_trials) * pairs_per_trial
  for (cond in conditions) {
    for (pid in seq_len(n_per_condition)) {
      pair_ids <- sample.int(noise_pool_size, draw_n,
                             replace = FALSE)
      trial_pairs <- matrix(pair_ids, nrow = n_trials,
                            ncol = pairs_per_trial)

      utilities_pos <- matrix(
        beta * proj[trial_pairs] +
          rgumbel(n_trials * pairs_per_trial),
        nrow = n_trials, ncol = pairs_per_trial
      )
      utilities_neg <- matrix(
        -beta * proj[trial_pairs] +
          rgumbel(n_trials * pairs_per_trial),
        nrow = n_trials, ncol = pairs_per_trial
      )

      utilities <- cbind(utilities_pos, utilities_neg)
      pick <- max.col(utilities, ties.method = "random")

      pair_idx_in_trial <- ifelse(pick <= pairs_per_trial,
                                  pick,
                                  pick - pairs_per_trial)
      response <- ifelse(pick <= pairs_per_trial, 1L, -1L)
      stim_id <- trial_pairs[
        cbind(seq_len(n_trials), pair_idx_in_trial)
      ]

      rts <- simulate_rts(
        n_trials,
        rt_contamination_fast,
        rt_contamination_slow
      )
      data_list[[k]] <- data.table::data.table(
        participant_id = sprintf("%s_%03d", cond, pid),
        condition      = cond,
        trial          = seq_len(n_trials),
        stimulus       = stim_id,
        response       = response,
        rt             = rts
      )
      k <- k + 1L
    }
  }
  data <- data.table::rbindlist(data_list)

  stimuli <- build_sim_stimuli(
    base_face  = base_face, p = pool$p, params = pool$params,
    img_size   = img_size, n_trials = noise_pool_size, seed = used_seed,
    noise_type = noise_type, nscales = nscales, sigma = sigma
  )
  rdata_path <- write_sim_rdata(
    stimuli = stimuli, dir = rdata_dir, method = "briefrc"
  )
  inform_rdata_written(rdata_path, persistent = !is.null(rdata_dir))

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  new_rcisignal_sim(
    data         = data,
    noise_matrix = pool$noise_matrix,
    base_face    = base_face,
    params       = pool$params,
    p            = pool$p,
    signal       = signal,
    rdata_path   = rdata_path,
    stimuli      = stimuli,
    meta         = list(
      method            = "briefrc",
      n_per_condition   = as.integer(n_per_condition),
      conditions        = as.character(conditions),
      n_trials          = as.integer(n_trials),
      images_per_trial  = as.integer(images_per_trial),
      noise_pool_size   = as.integer(noise_pool_size),
      img_size          = as.integer(img_size),
      signal_strength   = signal_strength,
      signal_region     = signal_region,
      seed              = used_seed,
      elapsed_secs      = elapsed
    )
  )
}

# ----- internals ------------------------------------------------------

#' @keywords internal
#' @noRd
ensure_rcicr <- function() {
  if (!requireNamespace("rcicr", quietly = TRUE)) {
    cli::cli_abort(c(
      "Simulating data requires the {.pkg rcicr} package for noise \\
       generation.",
      "i" = "Install: \\
             {.run remotes::install_github(\"rdotsch/rcicr\")}"
    ))
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
validate_simulate_args <- function(n_per_condition, conditions,
                                   n_trials, img_size,
                                   rt_contamination_fast,
                                   rt_contamination_slow) {
  if (!is.numeric(n_per_condition) || length(n_per_condition) != 1L ||
        n_per_condition < 1L) {
    cli::cli_abort("{.arg n_per_condition} must be a positive integer.")
  }
  if (!is.character(conditions) || length(conditions) < 1L ||
        anyDuplicated(conditions) > 0L) {
    cli::cli_abort(
      "{.arg conditions} must be a character vector of unique labels."
    )
  }
  if (!is.numeric(n_trials) || length(n_trials) != 1L ||
        n_trials < 2L) {
    cli::cli_abort("{.arg n_trials} must be an integer >= 2.")
  }
  if (!is.numeric(img_size) || length(img_size) != 1L ||
        img_size < 16L) {
    cli::cli_abort("{.arg img_size} must be an integer >= 16.")
  }
  for (nm in c("rt_contamination_fast", "rt_contamination_slow")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L ||
          val < 0 || val > 1) {
      cli::cli_abort("{.arg {nm}} must be in [0, 1].")
    }
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
resolve_signal_strength <- function(signal_strength) {
  if (is.numeric(signal_strength) && length(signal_strength) == 1L) {
    return(as.numeric(signal_strength))
  }
  if (is.character(signal_strength) &&
        length(signal_strength) == 1L) {
    return(switch(
      signal_strength,
      none   = 0,
      weak   = 0.5,
      strong = 2.0,
      cli::cli_abort(
        "{.arg signal_strength} string must be one of \\
         {.val none}, {.val weak}, {.val strong}."
      )
    ))
  }
  cli::cli_abort(
    "{.arg signal_strength} must be a single numeric or one of \\
     {.val none}, {.val weak}, {.val strong}."
  )
}

#' @keywords internal
#' @noRd
resolve_seed <- function(seed) {
  if (is.null(seed)) {
    return(sample.int(.Machine$integer.max, 1L))
  }
  if (!is.numeric(seed) || length(seed) != 1L) {
    cli::cli_abort("{.arg seed} must be a single integer or NULL.")
  }
  as.integer(seed)
}

#' @keywords internal
#' @noRd
save_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv())
  } else {
    NULL
  }
}

#' @keywords internal
#' @noRd
restore_rng_state <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(),
               inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
load_base_face <- function(base_image, img_size) {
  if (is.null(base_image)) {
    base_image <- system.file("extdata", "sim_base_face.png",
                              package = "rcisignal")
    if (!nzchar(base_image)) {
      cli::cli_abort(
        "Bundled base face not found at \\
         {.file inst/extdata/sim_base_face.png}."
      )
    }
  }
  if (is.matrix(base_image) && is.numeric(base_image)) {
    if (!all(dim(base_image) == c(img_size, img_size))) {
      cli::cli_abort(c(
        "{.arg base_image} matrix dimensions \\
         {.val {dim(base_image)}} do not match \\
         {.arg img_size} = {.val {img_size}}."
      ))
    }
    return(base_image)
  }
  if (!is.character(base_image) || length(base_image) != 1L) {
    cli::cli_abort(
      "{.arg base_image} must be NULL, a file path, or a numeric \\
       matrix."
    )
  }
  if (!file.exists(base_image)) {
    cli::cli_abort("Base image file not found: {.path {base_image}}")
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    cli::cli_abort(c(
      "Reading the base-image PNG requires the {.pkg png} package.",
      "i" = "Install: {.run install.packages(\"png\")}"
    ))
  }
  raw <- png::readPNG(base_image)
  if (length(dim(raw)) == 3L) {
    raw <- raw[, , 1L]
  }
  if (!all(dim(raw) == c(img_size, img_size))) {
    cli::cli_abort(c(
      "Base image at {.path {base_image}} is \\
       {.val {dim(raw)}} but {.arg img_size} = {.val {img_size}}.",
      "i" = "Provide a square PNG with side {.val {img_size}}, \\
             or set {.arg img_size} to {.val {dim(raw)[1L]}}."
    ))
  }
  raw
}

#' @keywords internal
#' @noRd
generate_noise_pool <- function(n_trials, img_size, noise_type,
                                nscales, sigma, progress) {
  p_basis <- rcicr::generateNoisePattern(
    img_size = img_size, nscales = nscales,
    noise_type = noise_type, sigma = sigma
  )
  ncoef <- guess_ncoef(p_basis)
  params <- matrix(stats::runif(n_trials * ncoef, -1, 1),
                   nrow = n_trials, ncol = ncoef)

  n_pix <- as.integer(img_size) * as.integer(img_size)
  noise_matrix <- matrix(NA_real_, nrow = n_pix, ncol = n_trials)

  if (isTRUE(progress)) {
    cli::cli_progress_bar(
      "Generating noise pool",
      total = n_trials,
      clear = TRUE
    )
  }
  for (i in seq_len(n_trials)) {
    noise_matrix[, i] <- as.vector(
      rcicr::generateNoiseImage(params[i, ], p_basis)
    )
    if (isTRUE(progress)) cli::cli_progress_update()
  }
  if (isTRUE(progress)) cli::cli_progress_done()

  list(
    noise_matrix = noise_matrix,
    params       = params,
    p            = p_basis
  )
}

#' @keywords internal
#' @noRd
guess_ncoef <- function(p_basis) {
  if (!is.null(p_basis$noise_type) &&
        identical(p_basis$noise_type, "gabor")) {
    return(length(p_basis$patches))
  }
  if (!is.null(p_basis$patchIdx)) {
    return(max(unlist(p_basis$patchIdx)))
  }
  if (!is.null(p_basis$patches)) {
    return(length(p_basis$patches))
  }
  cli::cli_abort(
    "Could not infer the number of basis coefficients from the \\
     {.fn rcicr::generateNoisePattern} return value."
  )
}

#' @keywords internal
#' @noRd
make_signal_vector <- function(img_size, signal_region) {
  mask <- make_face_mask(c(img_size, img_size), region = signal_region)
  as.numeric(mask)
}

#' @keywords internal
#' @noRd
simulate_rts <- function(n, frac_fast, frac_slow) {
  base_rt <- exp(stats::rnorm(n, mean = log(800), sd = 0.5)) + 150
  fast_n <- stats::rbinom(1L, n, frac_fast)
  slow_n <- stats::rbinom(1L, n, frac_slow)
  if (fast_n > 0L) {
    idx <- sample.int(n, fast_n)
    base_rt[idx] <- stats::runif(fast_n, 50, 200)
  }
  remaining <- setdiff(seq_len(n), if (fast_n > 0L) idx else integer())
  if (slow_n > 0L && length(remaining) > 0L) {
    slow_n <- min(slow_n, length(remaining))
    idx2 <- sample(remaining, slow_n)
    base_rt[idx2] <- stats::runif(slow_n, 5000, 20000)
  }
  round(base_rt)
}

#' Build the portable in-memory stimuli list
#'
#' Captures the rcicr-format stimuli env contents (base face,
#' basis, per-trial params, generator settings) as a plain list
#' that survives [saveRDS()]/[readRDS()]. Stored on the
#' `rcisignal_sim` object as `$stimuli` and accepted by every
#' rcicr-backed consumer (e.g. [ci_from_responses_2ifc()]) via a
#' `stimuli =` argument.
#'
#' @keywords internal
#' @noRd
build_sim_stimuli <- function(base_face, p, params, img_size,
                              n_trials, seed,
                              noise_type = "sinusoid",
                              nscales    = 5L,
                              sigma      = 25,
                              base_label = "base") {
  rcicr_ver <- tryCatch(
    as.character(utils::packageVersion("rcicr")),
    error = function(e) NA_character_
  )
  list(
    base_face         = base_face,
    base_label        = as.character(base_label),
    p                 = p,
    params            = params,
    img_size          = as.integer(img_size),
    n_trials          = as.integer(n_trials),
    seed              = as.integer(seed),
    noise_type        = as.character(noise_type),
    nscales           = as.integer(nscales),
    sigma             = as.numeric(sigma),
    generator_version = rcicr_ver
  )
}

#' Materialise a stimuli list to a fresh rcicr-format `.Rdata`
#'
#' Used by consumers (e.g. [ci_from_responses_2ifc()]) when the
#' user passes `stimuli =` instead of `rdata_path =`. Writes the
#' env contents to a fresh tempdir so calls into rcicr (which
#' `load()` a path) keep working unchanged.
#'
#' @keywords internal
#' @noRd
materialize_stimuli_rdata <- function(stimuli,
                                      method = c("2ifc",
                                                 "briefrc")) {
  method <- match.arg(method)
  required <- c("base_face", "base_label", "p", "params",
                "img_size", "n_trials", "seed", "noise_type",
                "nscales", "sigma")
  miss <- setdiff(required, names(stimuli))
  if (length(miss) > 0L) {
    cli::cli_abort(c(
      "{.arg stimuli} is missing required field{?s}: \\
       {.val {miss}}.",
      "i" = "Use the {.code $stimuli} element from a \\
             {.cls rcisignal_sim} object returned by \\
             {.fn simulate_2ifc_data} / \\
             {.fn simulate_briefrc_data}."
    ))
  }
  write_sim_rdata(stimuli = stimuli, dir = NULL, method = method)
}

#' Write an rcicr-format stimuli `.Rdata` to disk
#'
#' Single source of truth for the on-disk layout of the simulated
#' stimuli env. When `dir = NULL` writes to a fresh tempdir
#' (legacy behaviour). When `dir` is a character path the file
#' goes there under a stable filename (overwrites if present).
#'
#' @keywords internal
#' @noRd
write_sim_rdata <- function(stimuli, dir = NULL,
                            method = c("2ifc", "briefrc")) {
  method <- match.arg(method)
  base_label <- stimuli$base_label %||% "base"
  if (is.null(dir)) {
    dir <- tempfile("rcisignal_sim_")
  } else if (!is.character(dir) || length(dir) != 1L) {
    cli::cli_abort(
      "{.arg rdata_dir} must be a single string or NULL."
    )
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(dir, paste0(base_label, ".png"))
  if (requireNamespace("png", quietly = TRUE)) {
    png::writePNG(stimuli$base_face, png_path)
  } else {
    saveRDS(stimuli$base_face, sub("\\.png$", ".rds", png_path))
  }
  env <- new.env(parent = emptyenv())
  env$base_face_files     <- stats::setNames(list(png_path),
                                             base_label)
  env$base_faces          <- stats::setNames(list(stimuli$base_face),
                                             base_label)
  env$stimuli_params      <- stats::setNames(list(stimuli$params),
                                             base_label)
  env$img_size            <- as.integer(stimuli$img_size)
  env$n_trials            <- as.integer(stimuli$n_trials)
  env$seed                <- as.integer(stimuli$seed)
  env$label               <- "rcisignal_sim"
  env$stimulus_path       <- dir
  env$trial               <- seq_len(env$n_trials)
  env$generator_version   <- if (is.na(stimuli$generator_version %||%
                                         NA_character_)) {
    NA_character_
  } else {
    as.character(stimuli$generator_version)
  }
  env$use_same_parameters <- TRUE
  env$p                   <- stimuli$p
  env$noise_type          <- as.character(stimuli$noise_type)
  env$nscales             <- as.integer(stimuli$nscales)
  env$sigma               <- as.numeric(stimuli$sigma)
  fname <- paste0("rcisignal_sim_", method, "_stimuli.Rdata")
  rdata_file <- file.path(dir, fname)
  save(list = ls(env), envir = env, file = rdata_file)
  rdata_file
}

#' @keywords internal
#' @noRd
inform_rdata_written <- function(path, persistent) {
  if (persistent) {
    cli::cli_inform(c(
      "i" = "Wrote stimuli to {.path {path}}."
    ))
  } else {
    cli::cli_inform(c(
      "i" = "Wrote stimuli to {.path {path}} (session tempdir).",
      " " = "Pass {.arg rdata_dir} to persist across R sessions, \\
             or hand {.code $stimuli} to downstream consumers."
    ))
  }
  invisible(NULL)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Resolve `rdata` / `stimuli` to a usable path
#'
#' Consumers that delegate to rcicr need a file path (rcicr does
#' `load(rdata)` internally). When the user passes `stimuli`,
#' write a fresh tempdir-backed `.Rdata` and return its path so
#' the rest of the consumer pipeline is unchanged.
#'
#' @keywords internal
#' @noRd
resolve_rdata_input <- function(rdata_path, stimuli,
                                method = c("2ifc", "briefrc")) {
  method <- match.arg(method)
  if (!is.null(stimuli)) {
    if (!is.null(rdata_path)) {
      cli::cli_warn(
        "Both {.arg rdata_path} and {.arg stimuli} supplied; \\
         using {.arg stimuli}."
      )
    }
    return(materialize_stimuli_rdata(stimuli, method = method))
  }
  if (is.null(rdata_path)) {
    cli::cli_abort(c(
      "Pass either {.arg rdata_path} or {.arg stimuli}.",
      "i" = "{.arg stimuli} is the {.code $stimuli} element of an \\
             {.cls rcisignal_sim} object; use it when the \\
             {.code $rdata_path} on a saved sim no longer \\
             resolves (e.g. after {.fn saveRDS} / \\
             {.fn readRDS} across R sessions)."
    ))
  }
  if (!is.character(rdata_path) || length(rdata_path) != 1L) {
    cli::cli_abort(
      "{.arg rdata_path} must be a single string."
    )
  }
  if (!file.exists(rdata_path)) {
    cli::cli_abort(c(
      "rdata not found: {.path {rdata_path}}.",
      "i" = "If this path was stored on an {.cls rcisignal_sim} \\
             object before an R session restart, pass \\
             {.code stimuli = sim$stimuli} instead."
    ))
  }
  rdata_path
}

#' @keywords internal
#' @noRd
rgumbel <- function(n, location = 0, scale = 1) {
  u <- stats::runif(n, .Machine$double.eps, 1 - .Machine$double.eps)
  location - scale * log(-log(u))
}

