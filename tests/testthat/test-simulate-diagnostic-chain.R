## Regression test for the simulate_2ifc_data() -> downstream
## diagnostic chain. Catches recurrence of the v0.1.1 rdata-env
## bugs (CLAUDE.md §11.21): the rdata file written by
## simulate_2ifc_data() must satisfy both rcisignal's consumers
## (default baseimage = "base") and rcicr's reference-distribution
## regenerator (load(rdata); references noise_type).
##
## Gated on rcicr (Suggests) and skipped on CRAN-style runs.

skip_if_not_installed("rcicr")
testthat::skip_on_cran()

make_tiny_sim <- function(seed = 1L) {
  simulate_2ifc_data(
    n_per_condition = 5L,
    conditions      = c("target", "control"),
    n_trials        = 30L,
    img_size        = 256L,
    signal_strength = "weak",
    signal_region   = "eyes",
    seed            = seed,
    progress        = FALSE
  )
}

test_that("simulator + diagnose_infoval chain runs end-to-end", {
  sim <- make_tiny_sim()
  expect_true(file.exists(sim$rdata_path))

  res <- diagnose_infoval(
    sim$data,
    method   = "2ifc",
    rdata    = sim$rdata_path,
    iter     = 50L,
    seed     = 1L,
    progress = FALSE
  )
  expect_s3_class(res, "rcisignal_diag_result")
})

test_that("simulator + compute_infoval_summary chain runs end-to-end", {
  sim <- make_tiny_sim(seed = 2L)
  res <- compute_infoval_summary(
    sim$data,
    method = "2ifc",
    rdata  = sim$rdata_path,
    iter   = 50L
  )
  expect_s3_class(res, "rcisignal_diag_result")
})

test_that("simulator + check_rt_infoval_consistency chain runs", {
  sim <- make_tiny_sim(seed = 3L)
  res <- check_rt_infoval_consistency(
    sim$data,
    method = "2ifc",
    rdata  = sim$rdata_path,
    col_rt = "rt",
    iter   = 50L
  )
  expect_s3_class(res, "rcisignal_diag_result")
})

test_that("simulator + check_response_inversion chain runs", {
  sim <- make_tiny_sim(seed = 4L)
  res <- check_response_inversion(
    sim$data,
    method = "2ifc",
    rdata  = sim$rdata_path,
    iter   = 50L
  )
  expect_s3_class(res, "rcisignal_diag_result")
})

test_that("simulator + run_diagnostics chain runs end-to-end", {
  sim <- make_tiny_sim(seed = 5L)
  rep <- run_diagnostics(
    sim$data,
    method       = "2ifc",
    rdata        = sim$rdata_path,
    col_rt       = "rt",
    infoval_iter = 50L,
    seed         = 1L,
    progress     = FALSE
  )
  expect_s3_class(rep, "rcisignal_diag_report")
})
