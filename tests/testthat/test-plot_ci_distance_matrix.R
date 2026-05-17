## Tests for plot_ci_distance_matrix():
##   - basic run + return shape on group-mean vectors
##   - method = "raw" vs "normalised"
##   - auto-reduction of per-producer matrices
##   - mask + triangle + palette options
##   - file save (PNG + PDF)
##   - structural error paths

make_ci_vec <- function(seed, n_pix = 32L * 32L) {
  set.seed(seed)
  rnorm(n_pix) + sin(seq_len(n_pix) / 200)
}

n_side <- 32L
n_pix  <- n_side * n_side

test_that("plot_ci_distance_matrix runs and returns the right structure", {
  cis <- list(A = make_ci_vec(1), B = make_ci_vec(2),
              C = make_ci_vec(3))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_ci_distance_matrix(cis, img_dims = c(n_side, n_side))

  expect_type(out, "list")
  expect_named(out, c("distance_matrix", "distance_raw", "method",
                      "n_pixels_used", "mask", "palette", "triangle",
                      "file"))
  expect_true(is.matrix(out$distance_matrix))
  expect_equal(dim(out$distance_matrix), c(3L, 3L))
  expect_equal(rownames(out$distance_matrix), c("A", "B", "C"))
  expect_equal(diag(out$distance_matrix), c(A = 0, B = 0, C = 0))
  expect_equal(out$distance_matrix, t(out$distance_matrix))
  expect_equal(out$method, "raw")
  expect_equal(out$mask, "none")
  expect_equal(out$triangle, "full")
  expect_equal(out$palette, "viridis")
  expect_null(out$file)
  expect_equal(out$n_pixels_used, n_pix)
})

test_that("plot_ci_distance_matrix method = 'raw' makes the two matrices equal", {
  cis <- list(A = make_ci_vec(10), B = make_ci_vec(11),
              C = make_ci_vec(12))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_ci_distance_matrix(cis, img_dims = c(n_side, n_side),
                                 method = "raw")
  expect_equal(out$distance_matrix, out$distance_raw)
})

test_that("plot_ci_distance_matrix method = 'normalised' rescales by sqrt(n_pixels_used)", {
  cis <- list(A = make_ci_vec(20), B = make_ci_vec(21),
              C = make_ci_vec(22))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_ci_distance_matrix(cis, img_dims = c(n_side, n_side),
                                 mask   = "face",
                                 method = "normalised")
  expect_equal(out$distance_matrix,
               out$distance_raw / sqrt(out$n_pixels_used))
})

test_that("plot_ci_distance_matrix auto-reduces per-producer matrix to group mean", {
  v <- make_ci_vec(30)
  set.seed(31)
  per_producer <- matrix(rnorm(n_pix * 5L), n_pix, 5L)
  cis <- list(P = per_producer, V = v)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_ci_distance_matrix(cis, img_dims = c(n_side, n_side))
  expected <- sqrt(sum((rowMeans(per_producer) - v)^2))
  expect_equal(out$distance_raw["P", "V"], expected)
})

test_that("plot_ci_distance_matrix mask options decrease n_pixels_used", {
  cis <- list(A = make_ci_vec(40), B = make_ci_vec(41))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  none  <- plot_ci_distance_matrix(cis, img_dims = c(n_side, n_side),
                                   mask = "none")
  face  <- plot_ci_distance_matrix(cis, img_dims = c(n_side, n_side),
                                   mask = "face")
  upper <- plot_ci_distance_matrix(cis, img_dims = c(n_side, n_side),
                                   mask = "upper_face")
  lower <- plot_ci_distance_matrix(cis, img_dims = c(n_side, n_side),
                                   mask = "lower_face")

  expect_equal(none$n_pixels_used, n_pix)
  expect_lt(face$n_pixels_used,  none$n_pixels_used)
  expect_lt(upper$n_pixels_used, face$n_pixels_used)
  expect_lt(lower$n_pixels_used, face$n_pixels_used)
  expect_equal(upper$n_pixels_used + lower$n_pixels_used,
               face$n_pixels_used)
})

test_that("plot_ci_distance_matrix all four palettes run", {
  cis <- list(A = make_ci_vec(50), B = make_ci_vec(51))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in c("viridis", "inferno", "plasma", "rocket")) {
    expect_no_error(
      plot_ci_distance_matrix(cis, img_dims = c(n_side, n_side),
                              palette = p)
    )
  }
})

test_that("plot_ci_distance_matrix triangle = 'upper' / 'lower' both run", {
  cis <- list(A = make_ci_vec(60), B = make_ci_vec(61),
              C = make_ci_vec(62))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_ci_distance_matrix(cis,
                                          img_dims = c(n_side, n_side),
                                          triangle = "upper"))
  expect_no_error(plot_ci_distance_matrix(cis,
                                          img_dims = c(n_side, n_side),
                                          triangle = "lower"))
})

test_that("plot_ci_distance_matrix saves PNG via file =", {
  cis <- list(A = make_ci_vec(70), B = make_ci_vec(71),
              C = make_ci_vec(72))
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  out <- plot_ci_distance_matrix(cis, img_dims = c(n_side, n_side),
                                 file = tmp)
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 0)
  expect_equal(out$file, tmp)
})

test_that("plot_ci_distance_matrix saves PDF via file =", {
  cis <- list(A = make_ci_vec(80), B = make_ci_vec(81))
  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)
  out <- plot_ci_distance_matrix(cis, img_dims = c(n_side, n_side),
                                 file = tmp)
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 0)
})

test_that("plot_ci_distance_matrix rejects unsupported file extension", {
  cis <- list(A = make_ci_vec(90), B = make_ci_vec(91))
  tmp <- tempfile(fileext = ".jpg")
  expect_error(
    plot_ci_distance_matrix(cis, img_dims = c(n_side, n_side),
                            file = tmp),
    "\\.png.*\\.pdf"
  )
})

test_that("plot_ci_distance_matrix errors on bad input", {
  expect_error(plot_ci_distance_matrix("not a list"),
               "named list")
  expect_error(plot_ci_distance_matrix(list(A = rnorm(n_pix))),
               "at least")

  cis_mismatch <- list(A = rnorm(n_pix), B = rnorm(n_pix + 5L))
  expect_error(plot_ci_distance_matrix(cis_mismatch),
               "same number of pixels")
})
