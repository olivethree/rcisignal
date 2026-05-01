test_that("run_discriminability wraps cluster_test + dissimilarity in a report", {
  pair <- make_sig_pair(256L, 10L, seed = 1L)
  rep <- suppressWarnings(
    run_discriminability(pair$a, pair$b, img_dims = pair$img_dims,
                n_permutations = 50L, n_boot = 100L,
                seed = 1L, progress = FALSE)
  )
  expect_s3_class(rep, "rcisignal_rel_report")
  expect_setequal(names(rep$results),
                  c("cluster_test", "dissimilarity"))
  expect_s3_class(rep$results$cluster_test,  "rcisignal_rel_cluster_test")
  expect_s3_class(rep$results$dissimilarity, "rcisignal_rel_dissim")
  expect_equal(rep$method, "discriminability")
  expect_equal(rep$img_dims, pair$img_dims)
})
