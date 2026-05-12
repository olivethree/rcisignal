## These tests just call every S3 method to make sure they don't error.
## We don't assert on plot appearance.

test_that("print/summary/plot methods exist and do not error", {
  sig <- make_sig(64L, 10L)
  r <- suppressWarnings(
    rel_split_half(sig, n_permutations = 50L, progress = FALSE)
  )
  expect_invisible(print(r))
  expect_invisible(summary(r))
  expect_silent({ grDevices::pdf(NULL); plot(r); grDevices::dev.off() })

  lo <- suppressWarnings(rel_loo(sig))
  expect_invisible(print(lo))
  expect_invisible(summary(lo))
  expect_silent({ grDevices::pdf(NULL); plot(lo); grDevices::dev.off() })

  ic <- suppressWarnings(rel_icc(sig))
  expect_invisible(print(ic))
  expect_invisible(summary(ic))
  expect_silent({ grDevices::pdf(NULL); plot(ic); grDevices::dev.off() })

  pair <- make_sig_pair(64L, 8L)
  ct <- suppressWarnings(
    rel_cluster_test(pair$a, pair$b, img_dims = pair$img_dims,
                     n_permutations = 30L, progress = FALSE)
  )
  expect_invisible(print(ct))
  expect_invisible(summary(ct))
  expect_silent({ grDevices::pdf(NULL); plot(ct); grDevices::dev.off() })

  d <- suppressWarnings(
    rel_dissimilarity(pair$a, pair$b, n_boot = 50L, progress = FALSE)
  )
  expect_invisible(print(d))
  expect_invisible(summary(d))
  expect_silent({ grDevices::pdf(NULL); plot(d); grDevices::dev.off() })

  rep <- suppressWarnings(
    run_reliability(sig, n_permutations = 30L, progress = FALSE)
  )
  expect_invisible(print(rep))
  expect_invisible(summary(rep))
  expect_silent({ grDevices::pdf(NULL); plot(rep); grDevices::dev.off() })
})

test_that("cluster_test plot() accepts colour_bar = FALSE", {
  pair <- make_sig_pair(64L, 8L)
  ct <- suppressWarnings(
    rel_cluster_test(pair$a, pair$b, img_dims = pair$img_dims,
                     n_permutations = 30L, progress = FALSE)
  )
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(ct, colour_bar = FALSE))
  expect_no_error(plot(ct, colour_bar = TRUE))
})

test_that("plot.rcisignal_rel_pairwise_report draws a grid", {
  pa <- make_sig_pair(64L, 6L, seed = 1L)
  pb <- make_sig_pair(64L, 6L, seed = 2L)
  pc <- make_sig_pair(64L, 6L, seed = 3L)
  sigs <- list(A = pa$a, B = pa$b, C = pb$b)
  for (nm in names(sigs)) {
    attr(sigs[[nm]], "img_dims") <- pa$img_dims
  }
  pw <- suppressWarnings(
    run_discriminability_pairwise(
      sigs, n_permutations = 30L, n_boot = 50L,
      seed = 1L, progress = FALSE
    )
  )
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  mfrow_before <- graphics::par("mfrow")
  expect_no_error(plot(pw))
  expect_identical(graphics::par("mfrow"), mfrow_before)
})

test_that("plot.rcisignal_rel_pairwise_report warns when many pairs", {
  fake_child <- list(cluster_test = structure(
    list(
      observed_t        = rep(0, 16L),
      img_dims          = c(4L, 4L),
      method            = "threshold",
      cluster_threshold = 2.0,
      alpha             = 0.05,
      n_permutations    = 30L,
      pos_labels        = rep(0L, 16L),
      neg_labels        = rep(0L, 16L),
      clusters          = data.frame(
        cluster_id  = integer(0),
        direction   = character(0),
        significant = logical(0)
      )
    ),
    class = "rcisignal_rel_cluster_test"
  ))
  pairs <- paste0("p", seq_len(15L), "_vs_q", seq_len(15L))
  fake <- structure(
    list(
      results    = stats::setNames(replicate(15L, fake_child,
                                             simplify = FALSE), pairs),
      pairs      = data.frame(pair_id = pairs),
      conditions = letters[1:6],
      fwer       = "holm",
      alpha      = 0.05,
      method     = "threshold"
    ),
    class = "rcisignal_rel_pairwise_report"
  )
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_warning(plot(fake), "panels will be small")
  expect_no_warning(plot(fake, max_pairs = Inf))
})
