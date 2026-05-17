## Tests for plot_ci_correlogram():
##   - basic run + return shape on group-mean vectors
##   - auto-reduction of a per-producer signal matrix to a group mean
##   - triangle = "upper" / "lower" smoke tests
##   - mask = "face" / "upper_face" / "lower_face" smoke tests + n_pixels_used
##   - all three palette names run
##   - file = .png / .pdf save smoke tests
##   - extension validation
##   - structural error paths

make_ci_vec <- function(seed, n_pix = 32L * 32L) {
  set.seed(seed)
  rnorm(n_pix) + sin(seq_len(n_pix) / 200)
}

n_side <- 32L
n_pix  <- n_side * n_side

test_that("plot_ci_correlogram runs and returns the right structure", {
  cis <- list(A = make_ci_vec(1), B = make_ci_vec(2),
              C = make_ci_vec(3))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_ci_correlogram(cis, img_dims = c(n_side, n_side))
  expect_type(out, "list")
  expect_named(out, c("correlation_matrix", "n_pixels_used", "mask",
                      "palette", "triangle", "file"))
  expect_true(is.matrix(out$correlation_matrix))
  expect_equal(dim(out$correlation_matrix), c(3L, 3L))
  expect_equal(rownames(out$correlation_matrix), c("A", "B", "C"))
  expect_equal(out$correlation_matrix, t(out$correlation_matrix))
  expect_equal(diag(out$correlation_matrix), c(A = 1, B = 1, C = 1))
  expect_equal(out$mask, "none")
  expect_equal(out$triangle, "full")
  expect_equal(out$palette, "diverging")
  expect_null(out$file)
  expect_equal(out$n_pixels_used, n_pix)
})

test_that("plot_ci_correlogram auto-reduces per-producer matrix to group mean", {
  v <- make_ci_vec(10)
  set.seed(11)
  per_producer <- matrix(rnorm(n_pix * 5L), n_pix, 5L)
  cis <- list(P = per_producer, V = v)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_ci_correlogram(cis, img_dims = c(n_side, n_side))
  expected <- stats::cor(rowMeans(per_producer), v)
  expect_equal(out$correlation_matrix["P", "V"], expected)
})

test_that("plot_ci_correlogram triangle = 'upper' / 'lower' both run", {
  cis <- list(A = make_ci_vec(20), B = make_ci_vec(21),
              C = make_ci_vec(22))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_ci_correlogram(cis, img_dims = c(n_side, n_side),
                                      triangle = "upper"))
  expect_no_error(plot_ci_correlogram(cis, img_dims = c(n_side, n_side),
                                      triangle = "lower"))
})

test_that("plot_ci_correlogram mask options reduce n_pixels_used", {
  cis <- list(A = make_ci_vec(30), B = make_ci_vec(31))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  none  <- plot_ci_correlogram(cis, img_dims = c(n_side, n_side),
                               mask = "none")
  face  <- plot_ci_correlogram(cis, img_dims = c(n_side, n_side),
                               mask = "face")
  upper <- plot_ci_correlogram(cis, img_dims = c(n_side, n_side),
                               mask = "upper_face")
  lower <- plot_ci_correlogram(cis, img_dims = c(n_side, n_side),
                               mask = "lower_face")

  expect_equal(none$n_pixels_used, n_pix)
  expect_lt(face$n_pixels_used,  none$n_pixels_used)
  expect_lt(upper$n_pixels_used, face$n_pixels_used)
  expect_lt(lower$n_pixels_used, face$n_pixels_used)
  # upper + lower should partition the full-face mask
  expect_equal(upper$n_pixels_used + lower$n_pixels_used,
               face$n_pixels_used)
})

test_that("plot_ci_correlogram all three palettes run", {
  cis <- list(A = make_ci_vec(40), B = make_ci_vec(41))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_ci_correlogram(cis, img_dims = c(n_side, n_side),
                                      palette = "diverging"))
  expect_no_error(plot_ci_correlogram(cis, img_dims = c(n_side, n_side),
                                      palette = "diverging_puor"))
  expect_no_error(plot_ci_correlogram(cis, img_dims = c(n_side, n_side),
                                      palette = "diverging_brbg"))
})

test_that("plot_ci_correlogram saves PNG via file =", {
  cis <- list(A = make_ci_vec(50), B = make_ci_vec(51),
              C = make_ci_vec(52))
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  out <- plot_ci_correlogram(cis, img_dims = c(n_side, n_side),
                             file = tmp)
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 0)
  expect_equal(out$file, tmp)
})

test_that("plot_ci_correlogram saves PDF via file =", {
  cis <- list(A = make_ci_vec(60), B = make_ci_vec(61))
  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)
  out <- plot_ci_correlogram(cis, img_dims = c(n_side, n_side),
                             file = tmp)
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 0)
})

test_that("plot_ci_correlogram rejects unsupported file extension", {
  cis <- list(A = make_ci_vec(70), B = make_ci_vec(71))
  tmp <- tempfile(fileext = ".jpg")
  expect_error(
    plot_ci_correlogram(cis, img_dims = c(n_side, n_side), file = tmp),
    "\\.png.*\\.pdf"
  )
})

test_that("plot_ci_correlogram errors on bad input", {
  expect_error(plot_ci_correlogram("not a list"),
               "named list")
  expect_error(plot_ci_correlogram(list(A = rnorm(n_pix))),
               "at least two")
  expect_error(plot_ci_correlogram(list(rnorm(n_pix), rnorm(n_pix))),
               "fully named")

  cis_mismatch <- list(A = rnorm(n_pix), B = rnorm(n_pix + 5L))
  expect_error(plot_ci_correlogram(cis_mismatch),
               "same number of pixels")

  cis_bad_elem <- list(A = "hi", B = rnorm(n_pix))
  expect_error(plot_ci_correlogram(cis_bad_elem),
               "numeric vector or matrix")
})
