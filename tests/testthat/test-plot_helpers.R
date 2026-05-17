## Tests for the new plot helpers added in v0.2.x:
##   - plot_agreement_map()
##   - plot_dissimilarity_grid()
##
## These don't assert on appearance — only that the functions run
## without errors and return the expected structure.

test_that("plot_agreement_map runs and returns t_map of right length", {
  set.seed(1)
  n_side <- 32L
  n_pix  <- n_side * n_side
  rr <- row(matrix(0, n_side, n_side)) / n_side - 0.5
  cc <- col(matrix(0, n_side, n_side)) / n_side - 0.5
  d  <- sqrt(rr^2 + cc^2)
  mask_vec <- as.vector(pmax(0, 1 - d / 0.2))
  signal <- 0.5 * outer(mask_vec, runif(20L, 0.7, 1.3)) +
            matrix(rnorm(n_pix * 20L), n_pix, 20L)
  attr(signal, "img_dims") <- c(n_side, n_side)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- suppressWarnings(plot_agreement_map(signal))
  expect_type(out, "list")
  expect_named(out,
               c("t_map", "n", "img_dims", "mask", "zlim", "palette"))
  expect_length(out$t_map, n_pix)
  expect_equal(out$n, 20L)
  expect_equal(out$palette, "diverging")
  # Default diverging palette uses symmetric zlim around zero
  expect_equal(out$zlim[1L], -out$zlim[2L])
  # Signal pixels (centre) should have higher t than periphery
  centre_idx     <- which(mask_vec > 0.5)
  periphery_idx  <- which(mask_vec == 0)
  expect_gt(mean(out$t_map[centre_idx]),
            mean(out$t_map[periphery_idx]))
})

test_that("plot_agreement_map fire palette uses |t| and asymmetric zlim", {
  set.seed(3)
  n_side <- 32L
  n_pix  <- n_side * n_side
  signal <- matrix(rnorm(n_pix * 20L), n_pix, 20L)
  # Plant strong positive signal in one region, strong negative in another
  signal[1:200,         ] <- signal[1:200,         ] + 3
  signal[(n_pix - 199):n_pix, ] <- signal[(n_pix - 199):n_pix, ] - 3
  attr(signal, "img_dims") <- c(n_side, n_side)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- suppressWarnings(
    plot_agreement_map(signal, palette = "fire")
  )
  expect_equal(out$palette, "fire")
  # Asymmetric zlim starting at 0
  expect_equal(out$zlim[1L], 0)
  expect_gt(out$zlim[2L], 0)
  # The returned t_map is still SIGNED (display abs is only for plotting)
  expect_true(any(out$t_map > 0))
  expect_true(any(out$t_map < 0))
})

test_that("plot_agreement_map fire palette respects threshold and mask", {
  set.seed(4)
  signal <- matrix(rnorm(64L * 10L), 64L, 10L)
  attr(signal, "img_dims") <- c(8L, 8L)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  m <- rep(c(TRUE, FALSE), each = 32L)
  expect_no_error(
    plot_agreement_map(signal, palette = "fire",
                       mask = m, threshold = 0.5)
  )
})

test_that("plot_agreement_map rejects removed `viridis` palette option", {
  signal <- matrix(rnorm(64L * 10L), 64L, 10L)
  attr(signal, "img_dims") <- c(8L, 8L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(
    plot_agreement_map(signal, palette = "viridis"),
    "should be one of"
  )
})

test_that("plot_agreement_map respects mask and threshold args", {
  set.seed(2)
  signal <- matrix(rnorm(64L * 10L), 64L, 10L)
  attr(signal, "img_dims") <- c(8L, 8L)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  m <- rep(c(TRUE, FALSE), each = 32L)
  expect_no_error(
    plot_agreement_map(signal, mask = m, threshold = 0.5)
  )
})

test_that("plot_agreement_map errors on bad input", {
  expect_error(plot_agreement_map("not a matrix"),
               "numeric matrix")
  signal <- matrix(rnorm(50L * 5L), 50L, 5L)
  expect_error(plot_agreement_map(signal),
               "Cannot infer")  # not a perfect square, no img_dims
})

test_that("plot_dissimilarity_grid returns a tidy table and runs", {
  pair <- make_sig_pair(n_pix = 256L, n_p = 8L, seed = 1L)
  d1 <- suppressWarnings(rel_dissimilarity(pair$a, pair$b,
    n_boot = 100L, seed = 1L, progress = FALSE))
  d2 <- suppressWarnings(rel_dissimilarity(pair$a, pair$a,
    n_boot = 100L, seed = 2L, progress = FALSE))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_dissimilarity_grid("A vs B" = d1, "A vs A" = d2)
  expect_s3_class(out, "data.frame")
  expect_named(out,
               c("label", "observed", "ci_low", "ci_high"))
  expect_equal(nrow(out), 2L)
  # A-vs-A should have observed distance lower than A-vs-B
  expect_lt(out$observed[out$label == "A vs A"],
            out$observed[out$label == "A vs B"])
})

test_that("plot_agreement_map base_image = NULL keeps existing return shape", {
  set.seed(11)
  signal <- matrix(rnorm(64L * 10L), 64L, 10L)
  attr(signal, "img_dims") <- c(8L, 8L)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- suppressWarnings(plot_agreement_map(signal))
  expect_named(out,
               c("t_map", "n", "img_dims", "mask", "zlim", "palette"))
  expect_null(out$mask)
})

test_that("plot_agreement_map base_image as matrix runs and returns expected fields", {
  set.seed(12)
  signal <- matrix(rnorm(64L * 10L), 64L, 10L)
  attr(signal, "img_dims") <- c(8L, 8L)
  base <- matrix(0.5, 8L, 8L)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- suppressWarnings(
    plot_agreement_map(signal, base_image = base, threshold = 0.5)
  )
  expect_named(out,
               c("t_map", "n", "img_dims", "mask", "zlim", "palette"))
  expect_length(out$t_map, 64L)
  expect_equal(out$zlim[1L], -out$zlim[2L])

  # Same again on the fire palette to confirm both branches run.
  out2 <- suppressWarnings(
    plot_agreement_map(signal, base_image = base, palette = "fire",
                       threshold = 0.5)
  )
  expect_equal(out2$palette, "fire")
  expect_equal(out2$zlim[1L], 0)
})

test_that("plot_agreement_map base_image as PNG path runs", {
  skip_if_not_installed("png")
  png_path <- system.file("extdata", "sim_base_face.png",
                          package = "rcisignal")
  if (!nzchar(png_path) || !file.exists(png_path)) {
    skip("sim_base_face.png fixture not installed")
  }
  dims  <- dim(png::readPNG(png_path))[1:2]
  n_pix <- prod(dims)
  set.seed(13)
  signal <- matrix(rnorm(n_pix * 8L), n_pix, 8L)
  attr(signal, "img_dims") <- as.integer(dims)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(
    suppressWarnings(
      plot_agreement_map(signal, base_image = png_path)
    )
  )
})

test_that("plot_agreement_map errors when base_image dims disagree", {
  signal <- matrix(rnorm(64L * 10L), 64L, 10L)
  attr(signal, "img_dims") <- c(8L, 8L)
  bad_base <- matrix(0.5, 10L, 10L)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(
    suppressWarnings(
      plot_agreement_map(signal, base_image = bad_base)
    ),
    "do not match"
  )
})

test_that("plot_dissimilarity_grid errors on missing names or wrong class", {
  pair <- make_sig_pair(n_pix = 64L, n_p = 5L, seed = 1L)
  d <- suppressWarnings(rel_dissimilarity(pair$a, pair$b,
    n_boot = 50L, seed = 1L, progress = FALSE))
  expect_error(plot_dissimilarity_grid(d), "named")
  expect_error(plot_dissimilarity_grid("A" = list(foo = 1)),
               "rcisignal_rel_dissim")
  expect_error(plot_dissimilarity_grid(), "at least one")
})
