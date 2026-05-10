skip_if_no_rcicr <- function() {
  testthat::skip_if_not_installed("rcicr")
}

test_that("simulate_2ifc_data returns a well-formed rcisignal_sim", {
  skip_if_no_rcicr()
  sim <- simulate_2ifc_data(
    n_per_condition = 2L, n_trials = 8L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = "none", seed = 1L, progress = FALSE
  )
  expect_s3_class(sim, "rcisignal_sim")
  expect_s3_class(sim$data, "data.table")
  expect_named(sim$data,
               c("participant_id", "condition", "trial",
                 "stimulus", "response", "rt"))
  expect_equal(nrow(sim$data), 2L * 2L * 8L)
  expect_equal(dim(sim$noise_matrix), c(64L * 64L, 8L))
  expect_equal(dim(sim$base_face), c(64L, 64L))
  expect_true(all(sim$data$response %in% c(-1L, 1L)))
  expect_true(all(sim$data$rt > 0))
  expect_true(file.exists(sim$rdata_path))
  expect_equal(sim$meta$method, "2ifc")
  expect_equal(sim$meta$seed, 1L)
})

test_that("simulate_briefrc_data returns a well-formed rcisignal_sim", {
  skip_if_no_rcicr()
  sim <- simulate_briefrc_data(
    n_per_condition = 2L, n_trials = 5L, images_per_trial = 4L,
    img_size = 64L, base_image = matrix(0.5, 64, 64),
    signal_strength = "none", seed = 2L, progress = FALSE
  )
  expect_s3_class(sim, "rcisignal_sim")
  expect_equal(nrow(sim$data), 2L * 2L * 5L)
  expect_true(all(sim$data$response %in% c(-1L, 1L)))
  expect_equal(sim$meta$images_per_trial, 4L)
  expect_equal(sim$meta$noise_pool_size, 5L * 2L)
  expect_true(all(sim$data$stimulus %in%
                    seq_len(sim$meta$noise_pool_size)))
})

test_that("seeds are reproducible", {
  skip_if_no_rcicr()
  args <- list(
    n_per_condition = 2L, n_trials = 6L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = "weak", seed = 99L, progress = FALSE
  )
  a <- do.call(simulate_2ifc_data, args)
  b <- do.call(simulate_2ifc_data, args)
  expect_identical(a$data, b$data)
  expect_identical(a$noise_matrix, b$noise_matrix)
})

test_that("global RNG state is restored", {
  skip_if_no_rcicr()
  set.seed(123L)
  before <- .Random.seed
  simulate_2ifc_data(
    n_per_condition = 1L, n_trials = 4L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = "none", seed = 42L, progress = FALSE
  )
  expect_identical(.Random.seed, before)
})

test_that("signal_strength = 'none' yields ~50/50 responses (2IFC)", {
  skip_if_no_rcicr()
  sim <- simulate_2ifc_data(
    n_per_condition = 4L, n_trials = 200L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = "none", seed = 5L, progress = FALSE,
    rt_contamination_fast = 0, rt_contamination_slow = 0
  )
  prop_pos <- mean(sim$data$response == 1L)
  expect_gt(prop_pos, 0.40)
  expect_lt(prop_pos, 0.60)
})

test_that("strong signal flips response distribution", {
  skip_if_no_rcicr()
  sim_none <- simulate_2ifc_data(
    n_per_condition = 4L, n_trials = 100L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = "none", seed = 11L, progress = FALSE
  )
  sim_strong <- simulate_2ifc_data(
    n_per_condition = 4L, n_trials = 100L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = 8, seed = 11L, progress = FALSE
  )
  agree_strong <- mean(
    sim_strong$data[, list(p_pos = mean(response == 1L)),
                    by = "participant_id"]$p_pos
  )
  agree_none <- mean(
    sim_none$data[, list(p_pos = mean(response == 1L)),
                  by = "participant_id"]$p_pos
  )
  expect_gt(abs(agree_strong - 0.5), abs(agree_none - 0.5))
})

test_that("invalid args abort", {
  skip_if_no_rcicr()
  expect_error(
    simulate_2ifc_data(n_per_condition = 0L, progress = FALSE),
    "n_per_condition"
  )
  expect_error(
    simulate_2ifc_data(conditions = c("a", "a"), progress = FALSE),
    "conditions"
  )
  expect_error(
    simulate_2ifc_data(n_trials = 1L, progress = FALSE),
    "n_trials"
  )
  expect_error(
    simulate_2ifc_data(rt_contamination_fast = -0.1,
                       progress = FALSE),
    "rt_contamination_fast"
  )
  expect_error(
    simulate_2ifc_data(signal_strength = "bogus",
                       progress = FALSE),
    "signal_strength"
  )
  expect_error(
    simulate_briefrc_data(images_per_trial = 3L, progress = FALSE),
    "images_per_trial"
  )
})

test_that("simulated 2IFC data flows through run_diagnostics", {
  skip_if_no_rcicr()
  sim <- simulate_2ifc_data(
    n_per_condition = 2L, n_trials = 10L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = "weak", seed = 17L, progress = FALSE
  )
  rep <- run_diagnostics(sim$data, method = "2ifc", col_rt = "rt")
  expect_s3_class(rep, "rcisignal_diag_report")
})

test_that("simulated Brief-RC data flows through run_diagnostics", {
  skip_if_no_rcicr()
  sim <- simulate_briefrc_data(
    n_per_condition = 2L, n_trials = 8L, images_per_trial = 4L,
    img_size = 64L, base_image = matrix(0.5, 64, 64),
    signal_strength = "weak", seed = 18L, progress = FALSE
  )
  rep <- run_diagnostics(
    sim$data, method = "briefrc",
    noise_matrix = sim$noise_matrix, col_rt = "rt"
  )
  expect_s3_class(rep, "rcisignal_diag_report")
})

test_that("rdata_dir persists the stimuli file at a stable path", {
  skip_if_no_rcicr()
  tmp <- withr::local_tempdir()
  sim <- suppressMessages(simulate_2ifc_data(
    n_per_condition = 1L, n_trials = 6L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = "none", rdata_dir = tmp,
    seed = 1L, progress = FALSE
  ))
  expect_true(file.exists(sim$rdata_path))
  expect_equal(normalizePath(dirname(sim$rdata_path)),
               normalizePath(tmp))
  expect_match(basename(sim$rdata_path),
               "^rcisignal_sim_2ifc_stimuli\\.Rdata$")

  sim_br <- suppressMessages(simulate_briefrc_data(
    n_per_condition = 1L, n_trials = 4L, images_per_trial = 4L,
    img_size = 64L, base_image = matrix(0.5, 64, 64),
    signal_strength = "none", rdata_dir = tmp,
    seed = 2L, progress = FALSE
  ))
  expect_true(file.exists(sim_br$rdata_path))
  expect_match(basename(sim_br$rdata_path),
               "^rcisignal_sim_briefrc_stimuli\\.Rdata$")
})

test_that("rdata_dir = NULL keeps legacy tempdir behaviour", {
  skip_if_no_rcicr()
  sim <- suppressMessages(simulate_2ifc_data(
    n_per_condition = 1L, n_trials = 6L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = "none", seed = 1L, progress = FALSE
  ))
  expect_true(file.exists(sim$rdata_path))
  expect_true(grepl(normalizePath(tempdir(), mustWork = FALSE),
                    normalizePath(sim$rdata_path),
                    fixed = TRUE))
})

test_that("simulator emits an informational write message", {
  skip_if_no_rcicr()
  expect_message(
    simulate_2ifc_data(
      n_per_condition = 1L, n_trials = 4L, img_size = 64L,
      base_image = matrix(0.5, 64, 64),
      signal_strength = "none", seed = 1L, progress = FALSE
    ),
    "Wrote stimuli to"
  )
})

test_that("sim$stimuli is self-contained and round-trips via saveRDS", {
  skip_if_no_rcicr()
  sim <- suppressMessages(simulate_2ifc_data(
    n_per_condition = 1L, n_trials = 6L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = "none", seed = 1L, progress = FALSE
  ))
  expect_type(sim$stimuli, "list")
  expect_true(all(c("base_face", "p", "params", "img_size", "n_trials",
                    "seed", "noise_type", "nscales", "sigma",
                    "base_label") %in% names(sim$stimuli)))

  rds <- tempfile(fileext = ".rds")
  saveRDS(sim, rds)
  file.remove(sim$rdata_path)
  expect_false(file.exists(sim$rdata_path))

  restored <- readRDS(rds)
  expect_type(restored$stimuli, "list")
  expect_equal(dim(restored$stimuli$base_face), c(64L, 64L))
})

test_that("ci_from_responses_2ifc accepts stimuli= without rdata_path", {
  skip_if_not_installed("rcicr")
  skip_if_not_installed("foreach")
  skip_if_not_installed("tibble")
  skip_if_not_installed("dplyr")
  sim <- suppressMessages(simulate_2ifc_data(
    n_per_condition = 2L, n_trials = 8L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = "none", seed = 1L, progress = FALSE
  ))
  # Simulate the path-vanishing scenario.
  unlink(dirname(sim$rdata_path), recursive = TRUE)
  expect_false(file.exists(sim$rdata_path))

  res <- ci_from_responses_2ifc(sim$data, stimuli = sim$stimuli)
  expect_equal(dim(res$signal_matrix), c(64L * 64L, 2L * 2L))
})

test_that("simulate_2ifc_data writes a base-face PNG and exposes its path", {
  skip_if_no_rcicr()
  skip_if_not_installed("png")
  tmp <- withr::local_tempdir()
  sim <- suppressMessages(simulate_2ifc_data(
    n_per_condition = 1L, n_trials = 4L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = "none", rdata_dir = tmp,
    seed = 1L, progress = FALSE
  ))
  expect_true(file.exists(sim$base_image_path))
  expect_equal(normalizePath(dirname(sim$base_image_path)),
               normalizePath(tmp))
  expect_match(basename(sim$base_image_path),
               "^rcisignal_sim_2ifc_base_face\\.png$")
  img <- png::readPNG(sim$base_image_path)
  if (length(dim(img)) == 3L) img <- img[, , 1L]
  expect_equal(dim(img), dim(sim$base_face))
})

test_that("simulate_briefrc_data writes a base-face PNG and exposes its path", {
  skip_if_no_rcicr()
  skip_if_not_installed("png")
  tmp <- withr::local_tempdir()
  sim <- suppressMessages(simulate_briefrc_data(
    n_per_condition = 1L, n_trials = 4L, images_per_trial = 4L,
    img_size = 64L, base_image = matrix(0.5, 64, 64),
    signal_strength = "none", rdata_dir = tmp,
    seed = 2L, progress = FALSE
  ))
  expect_true(file.exists(sim$base_image_path))
  expect_equal(normalizePath(dirname(sim$base_image_path)),
               normalizePath(tmp))
  expect_match(basename(sim$base_image_path),
               "^rcisignal_sim_briefrc_base_face\\.png$")
  img <- png::readPNG(sim$base_image_path)
  if (length(dim(img)) == 3L) img <- img[, , 1L]
  expect_equal(dim(img), dim(sim$base_face))
})

test_that("simulator base-face PNG also written under tempdir mode", {
  skip_if_no_rcicr()
  skip_if_not_installed("png")
  sim <- suppressMessages(simulate_2ifc_data(
    n_per_condition = 1L, n_trials = 4L, img_size = 64L,
    base_image = matrix(0.5, 64, 64),
    signal_strength = "none", seed = 1L, progress = FALSE
  ))
  expect_true(file.exists(sim$base_image_path))
  expect_equal(normalizePath(dirname(sim$base_image_path)),
               normalizePath(dirname(sim$rdata_path)))
})

test_that("resolve_rdata_input errors when both inputs are missing", {
  skip_if_not_installed("rcicr")
  resp <- data.frame(
    participant_id = c("p1", "p1"),
    stimulus       = c(1L, 2L),
    response       = c(-1L, 1L)
  )
  expect_error(
    ci_from_responses_2ifc(resp),
    "rdata_path"
  )
})
