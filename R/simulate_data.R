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
#' @param seed Integer or `NULL`. When `NULL`, a random seed is
#'   drawn and stored on the result so the run is reproducible.
#' @param progress Logical. Show a `cli` progress bar during noise
#'   generation. Default `TRUE`.
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
#'   * `rdata_path` — path to an rcicr-format `.Rdata` file written
#'     to a session tempdir, suitable for
#'     [ci_from_responses_2ifc()] / [compute_infoval_summary()]
#'     and other downstream functions that take an `rdata` argument.
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
#' sim <- simulate_2ifc_data(n_per_condition = 5, n_trials = 50)
#' report <- run_diagnostics(sim$data, method = "2ifc", col_rt = "rt")
#' print(report)
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

  rdata_path <- write_sim_rdata_2ifc(
    base_face = base_face, p = pool$p, params = pool$params,
    img_size = img_size, n_trials = n_trials, seed = used_seed
  )

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  new_rcisignal_sim(
    data         = data,
    noise_matrix = pool$noise_matrix,
    base_face    = base_face,
    params       = pool$params,
    p            = pool$p,
    signal       = signal,
    rdata_path   = rdata_path,
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
#' @param n_per_condition,conditions,n_trials,img_size,base_image,signal_strength,signal_region,rt_contamination_fast,rt_contamination_slow,noise_type,nscales,sigma,seed,progress
#'   See [simulate_2ifc_data()]. Note: `n_trials` here means the
#'   number of Brief-RC trials per participant, not the noise pool
#'   size. Default `n_trials = 500`.
#' @param images_per_trial Integer (even). Number of images shown
#'   per trial; half are original and half are inverted versions of
#'   the same noise patterns. Default `12` (= 6 pairs).
#' @param noise_pool_size Integer. Total number of noise patterns to
#'   pre-generate. Default `n_trials * (images_per_trial / 2)`,
#'   i.e. enough so each participant samples without replacement.
#'   Pass a larger value to study sub-sampling.
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
#' sim <- simulate_briefrc_data(n_per_condition = 5, n_trials = 30)
#' report <- run_diagnostics(sim$data, method = "briefrc",
#'                           noise_matrix = sim$noise_matrix,
#'                           col_rt = "rt")
#' print(report)
#' }
#' @export
simulate_briefrc_data <- function(n_per_condition       = 50L,
                                  conditions            = c("target",
                                                            "control"),
                                  n_trials              = 500L,
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
                                  seed                  = NULL,
                                  progress              = TRUE) {
  ensure_rcicr()
  validate_simulate_args(
    n_per_condition, conditions, n_trials, img_size,
    rt_contamination_fast, rt_contamination_slow
  )
  if (!is.numeric(images_per_trial) || length(images_per_trial) != 1L ||
        images_per_trial < 2L || images_per_trial %% 2L != 0L) {
    cli::cli_abort(
      "{.arg images_per_trial} must be an even integer >= 2."
    )
  }
  pairs_per_trial <- as.integer(images_per_trial) %/% 2L
  if (is.null(noise_pool_size)) {
    noise_pool_size <- as.integer(n_trials) * pairs_per_trial
  }
  if (!is.numeric(noise_pool_size) || length(noise_pool_size) != 1L ||
        noise_pool_size < pairs_per_trial) {
    cli::cli_abort(c(
      "{.arg noise_pool_size} must be a single integer \\
       >= {.val {pairs_per_trial}}.",
      "i" = "Pool must contain at least one trial's worth of pairs."
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

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  new_rcisignal_sim(
    data         = data,
    noise_matrix = pool$noise_matrix,
    base_face    = base_face,
    params       = pool$params,
    p            = pool$p,
    signal       = signal,
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

#' @keywords internal
#' @noRd
write_sim_rdata_2ifc <- function(base_face, p, params, img_size,
                                 n_trials, seed,
                                 base_label = "sim_base") {
  dir <- tempfile("rcisignal_sim_")
  dir.create(dir, recursive = TRUE)
  png_path <- file.path(dir, paste0(base_label, ".png"))
  if (requireNamespace("png", quietly = TRUE)) {
    png::writePNG(base_face, png_path)
  } else {
    saveRDS(base_face, sub("\\.png$", ".rds", png_path))
  }
  env <- new.env(parent = emptyenv())
  env$base_face_files     <- stats::setNames(list(png_path),
                                             base_label)
  env$base_faces          <- stats::setNames(list(base_face),
                                             base_label)
  env$stimuli_params      <- stats::setNames(list(params),
                                             base_label)
  env$img_size            <- as.integer(img_size)
  env$n_trials            <- as.integer(n_trials)
  env$seed                <- as.integer(seed)
  env$label               <- "rcisignal_sim"
  env$stimulus_path       <- dir
  env$trial               <- seq_len(env$n_trials)
  env$generator_version   <- as.character(
    utils::packageVersion("rcicr")
  )
  env$use_same_parameters <- TRUE
  env$p                   <- p
  rdata_file <- file.path(dir, "rcisignal_sim_stimuli.Rdata")
  save(list = ls(env), envir = env, file = rdata_file)
  rdata_file
}

#' @keywords internal
#' @noRd
rgumbel <- function(n, location = 0, scale = 1) {
  u <- stats::runif(n, .Machine$double.eps, 1 - .Machine$double.eps)
  location - scale * log(-log(u))
}

