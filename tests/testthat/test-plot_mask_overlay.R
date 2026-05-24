## Smoke + input-handling tests for plot_mask_overlay().
## These don't assert on appearance — only that the function runs,
## returns invisibly, restores par() on exit, and handles each of
## the accepted (base, mask) input combinations.

open_null_device <- function() {
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
}

test_that("plot_mask_overlay accepts matrix base + vector mask", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  base <- matrix(stats::runif(64L * 64L), 64L, 64L)
  mask <- make_face_mask(c(64L, 64L), region = "eyes")
  res <- plot_mask_overlay(base, mask)
  expect_null(res)
})

test_that("plot_mask_overlay accepts matrix base + matrix mask", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  base <- matrix(stats::runif(32L * 32L), 32L, 32L)
  mask_mat <- matrix(make_face_mask(c(32L, 32L), region = "mouth"),
                     nrow = 32L, ncol = 32L)
  expect_no_error(plot_mask_overlay(base, mask_mat))
})

test_that("plot_mask_overlay accepts a path base + path mask", {
  skip_if_not_installed("png")
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  fixture <- system.file("extdata", "sim_base_face.png",
                         package = "rcisignal")
  skip_if(!nzchar(fixture) || !file.exists(fixture),
          "sim_base_face.png fixture unavailable")

  dims <- dim(png::readPNG(fixture))[1:2]
  mask_path <- tempfile(fileext = ".png")
  oval <- make_face_mask(dims, region = "full")
  png::writePNG(matrix(as.numeric(oval), dims[1L], dims[2L]),
                mask_path)

  expect_no_error(plot_mask_overlay(fixture, mask_path))
})

test_that("plot_mask_overlay aborts on dimension mismatch", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  base <- matrix(0.5, 32L, 32L)
  expect_error(
    plot_mask_overlay(base, rep(TRUE, 100L)),
    "length does not match"
  )
  bad_mat <- matrix(TRUE, 16L, 16L)
  expect_error(
    plot_mask_overlay(base, bad_mat),
    "dimensions do not match"
  )
})

test_that("plot_mask_overlay aborts on bad alpha / overlay_col / base", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  base <- matrix(0.5, 16L, 16L)
  mask <- rep(c(TRUE, FALSE), 128L)
  expect_error(plot_mask_overlay(base, mask, alpha = 2),
               "in \\[0, 1\\]")
  expect_error(plot_mask_overlay(base, mask, overlay_col = c("red", "blue")),
               "single color")
  expect_error(plot_mask_overlay("does/not/exist.png", mask),
               "not found")
  expect_error(plot_mask_overlay(list(1, 2), mask),
               "numeric matrix or a path")
})

test_that("plot_mask_overlay accepts a region shortcut", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  base <- matrix(0.5, 48L, 48L)
  expect_no_error(plot_mask_overlay(base, region = "left_eye"))
  expect_no_error(plot_mask_overlay(base, region = "eyes",
                                    region_bounds = c(0.30, 0.50, 0.10, 0.90)))

  expect_error(plot_mask_overlay(base), "Supply one of")
  expect_error(
    plot_mask_overlay(base, mask = rep(TRUE, 48L * 48L),
                      region = "eyes"),
    "either"
  )
  expect_error(
    plot_mask_overlay(base, mask = rep(TRUE, 48L * 48L),
                      region_bounds = c(0.30, 0.50, 0.10, 0.90)),
    "only valid with"
  )
})

test_that("plot_mask_overlay restores par() on exit", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  base <- matrix(0.5, 24L, 24L)
  mask <- make_face_mask(c(24L, 24L), region = "nose")

  mfrow_before <- graphics::par("mfrow")
  pty_before   <- graphics::par("pty")
  plot_mask_overlay(base, mask)
  expect_identical(graphics::par("mfrow"), mfrow_before)
  expect_identical(graphics::par("pty"),   pty_before)
})
