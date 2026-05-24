## Verify every stage-1 function rejects a stage-2 (group-averaged)
## object with the teaching message.

make_group_ci <- function(n_pix = 64L) {
  sm <- make_sig(n_pix = n_pix, n_p = 4L, seed = 1L)
  colnames(sm) <- sprintf("p%02d", 1:4)
  responses <- data.frame(
    participant_id = sprintf("p%02d", 1:4),
    condition      = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  group_ci(sm, responses, by = "condition")
}

test_that("infoval rejects rcisignal_group_ci", {
  gcis <- make_group_ci()
  noise <- matrix(rnorm(64L * 50L), 64L, 50L)
  expect_error(
    infoval(gcis, noise, trial_counts = c(A = 30L, B = 30L),
            iter = 50L, progress = FALSE),
    regexp = "stage 1"
  )
})

test_that("rel_split_half rejects rcisignal_group_ci", {
  expect_error(
    rel_split_half(make_group_ci(), n_permutations = 50L,
                   progress = FALSE),
    regexp = "stage 1"
  )
})

test_that("rel_icc rejects rcisignal_group_ci", {
  expect_error(
    rel_icc(make_group_ci()),
    regexp = "stage 1"
  )
})

test_that("rel_loo rejects rcisignal_group_ci", {
  expect_error(
    rel_loo(make_group_ci()),
    regexp = "stage 1"
  )
})

test_that("rel_cluster_test rejects rcisignal_group_ci (either arg)", {
  ok  <- make_sig(n_pix = 64L, n_p = 4L, seed = 2L)
  bad <- make_group_ci()
  expect_error(
    rel_cluster_test(bad, ok, img_dims = c(8L, 8L),
                     n_permutations = 50L, progress = FALSE),
    regexp = "stage 1"
  )
  expect_error(
    rel_cluster_test(ok, bad, img_dims = c(8L, 8L),
                     n_permutations = 50L, progress = FALSE),
    regexp = "stage 1"
  )
})

test_that("rel_dissimilarity rejects rcisignal_group_ci (either arg)", {
  ok  <- make_sig(n_pix = 64L, n_p = 4L, seed = 3L)
  bad <- make_group_ci()
  expect_error(
    rel_dissimilarity(bad, ok, n_boot = 50L, progress = FALSE),
    regexp = "stage 1"
  )
  expect_error(
    rel_dissimilarity(ok, bad, n_boot = 50L, progress = FALSE),
    regexp = "stage 1"
  )
})

test_that("agreement_map_test rejects rcisignal_group_ci", {
  expect_error(
    agreement_map_test(make_group_ci(), n_permutations = 50L,
                       progress = FALSE),
    regexp = "stage 1"
  )
})

test_that("pixel_t_test rejects rcisignal_group_ci (either arg)", {
  ok  <- make_sig(n_pix = 64L, n_p = 4L, seed = 4L)
  bad <- make_group_ci()
  expect_error(pixel_t_test(bad, ok), regexp = "stage 1")
  expect_error(pixel_t_test(ok, bad), regexp = "stage 1")
})

test_that("run_reliability rejects rcisignal_group_ci", {
  expect_error(
    run_reliability(make_group_ci(), n_permutations = 50L,
                    progress = FALSE),
    regexp = "stage 1"
  )
})

test_that("run_discriminability rejects rcisignal_group_ci", {
  ok  <- make_sig(n_pix = 64L, n_p = 4L, seed = 5L)
  bad <- make_group_ci()
  expect_error(
    run_discriminability(bad, ok, img_dims = c(8L, 8L),
                         n_permutations = 50L, n_boot = 50L,
                         progress = FALSE),
    regexp = "stage 1"
  )
  expect_error(
    run_discriminability(ok, bad, img_dims = c(8L, 8L),
                         n_permutations = 50L, n_boot = 50L,
                         progress = FALSE),
    regexp = "stage 1"
  )
})

test_that("run_discriminability_pairwise rejects rcisignal_group_ci", {
  ok  <- make_sig(n_pix = 64L, n_p = 4L, seed = 6L)
  bad <- make_group_ci()
  expect_error(
    run_discriminability_pairwise(
      list(A = ok, B = bad),
      img_dims = c(8L, 8L),
      n_permutations = 50L, n_boot = 50L, progress = FALSE
    ),
    regexp = "stage 1"
  )
})
