test_that("face_mask returns a logical vector of the right length", {
  m <- make_face_mask(c(64L, 64L))
  expect_type(m, "logical")
  expect_length(m, 64L * 64L)
})

test_that("face_mask 'full' covers ~half the image at default geometry", {
  m <- make_face_mask(c(128L, 128L))
  cov <- mean(m)
  # Default oval has hw=0.35, hh=0.45, area = pi*hw*hh ~= 0.495
  expect_gt(cov, 0.40)
  expect_lt(cov, 0.55)
})

test_that("elliptical sub-regions are subsets of the full oval", {
  full <- make_face_mask(c(128L, 128L), region = "full")
  for (region in c("nose", "mouth", "upper_face", "lower_face")) {
    sub <- make_face_mask(c(128L, 128L), region = region)
    expect_true(all(sub <= full),
                info = sprintf("region=%s leaks outside full oval", region))
    expect_true(any(sub),
                info = sprintf("region=%s is empty", region))
  }
})

test_that("rectangle eye regions are non-empty and independent of oval", {
  for (region in c("eyes", "left_eye", "right_eye")) {
    sub <- make_face_mask(c(128L, 128L), region = region)
    expect_true(any(sub),
                info = sprintf("region=%s is empty", region))
  }
  # eyes covers (and union of left+right is a subset of) the wide band
  eyes  <- make_face_mask(c(128L, 128L), region = "eyes")
  left  <- make_face_mask(c(128L, 128L), region = "left_eye")
  right <- make_face_mask(c(128L, 128L), region = "right_eye")
  expect_true(all(left  <= eyes))
  expect_true(all(right <= eyes))
  expect_false(any(left & right))
})

test_that("region_bounds tunes rectangle regions and rejects misuse", {
  m1 <- make_face_mask(c(64L, 64L), region = "left_eye")
  m2 <- make_face_mask(c(64L, 64L), region = "left_eye",
                       region_bounds = c(0.40, 0.50, 0.30, 0.50))
  expect_false(identical(m1, m2))
  expect_true(any(m2))

  # rejects bounds on non-rectangle regions
  expect_error(
    make_face_mask(c(64L, 64L), region = "mouth",
                   region_bounds = c(0.40, 0.50, 0.30, 0.50)),
    "rectangle"
  )

  # validates shape
  expect_error(
    make_face_mask(c(64L, 64L), region = "eyes",
                   region_bounds = c(0.40, 0.50, 0.30)),
    "length 4"
  )
  # rejects out-of-range
  expect_error(
    make_face_mask(c(64L, 64L), region = "eyes",
                   region_bounds = c(-0.10, 0.50, 0.30, 0.50)),
    "\\[0, 1\\]"
  )
  # rejects unordered pairs
  expect_error(
    make_face_mask(c(64L, 64L), region = "eyes",
                   region_bounds = c(0.50, 0.40, 0.30, 0.50)),
    "row_min < row_max"
  )
})

test_that("face_mask accepts a single integer for square images", {
  expect_identical(
    make_face_mask(64L),
    make_face_mask(c(64L, 64L))
  )
})

test_that("face_mask rejects bad img_dims", {
  expect_error(make_face_mask(c(0L, 64L)), "positive")
  expect_error(make_face_mask(c(64L, 64L, 64L)), "positive")
})

test_that("plot_face_mask accepts a region shortcut", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(out <- plot_face_mask(region = "left_eye",
                                      img_dims = c(48L, 48L)))
  expect_true(is.matrix(out))
  expect_true(any(out))
  expect_error(plot_face_mask(), "Supply one of")
  expect_error(
    plot_face_mask(region = "left_eye"),
    "img_dims"
  )
  m <- make_face_mask(c(16L, 16L), region = "full")
  expect_error(
    plot_face_mask(m, region = "left_eye", img_dims = c(16L, 16L)),
    "either"
  )
})

test_that("plot_face_mask renders for vector + img_dims", {
  m <- make_face_mask(c(64L, 64L), region = "eyes")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(out <- plot_face_mask(m, img_dims = c(64L, 64L)))
  expect_true(is.matrix(out))
  expect_true(is.logical(out))
  expect_identical(dim(out), c(64L, 64L))
  expect_equal(sum(out), sum(m))
})

test_that("plot_face_mask accepts a logical matrix directly", {
  mat <- matrix(make_face_mask(c(48L, 48L)), 48L, 48L)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot_face_mask(mat))
})

test_that("plot_face_mask requires img_dims for vector input", {
  m <- make_face_mask(c(32L, 32L))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_error(plot_face_mask(m), "img_dims")
})

test_that("plot_face_mask warns and drops mismatched base_image", {
  m <- make_face_mask(c(32L, 32L))
  bad_base <- matrix(0.5, 16L, 16L)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_warning(
    plot_face_mask(m, img_dims = c(32L, 32L), base_image = bad_base),
    "do not match"
  )
})

test_that("plot_face_mask validates alpha", {
  m <- make_face_mask(c(32L, 32L))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_error(plot_face_mask(m, img_dims = c(32L, 32L), alpha = 2),
               "alpha")
})
