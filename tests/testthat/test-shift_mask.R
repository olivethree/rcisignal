test_that("shift_mask preserves count for an in-bounds shift", {
  m <- make_face_mask(c(64L, 64L), region = "mouth")
  shifted <- shift_mask(m, vertical = 4L, img_dims = c(64L, 64L))
  expect_type(shifted, "logical")
  expect_length(shifted, 64L * 64L)
  expect_equal(sum(m), sum(shifted))
})

test_that("shift_mask matrix-in returns matrix-out of same dim", {
  mat <- matrix(make_face_mask(c(48L, 48L), region = "nose"),
                48L, 48L)
  out <- shift_mask(mat, vertical = 2L, horizontal = 3L)
  expect_true(is.matrix(out))
  expect_true(is.logical(out))
  expect_identical(dim(out), dim(mat))
})

test_that("shift_mask drops pixels that fall off the image", {
  mat <- matrix(FALSE, 16L, 16L)
  mat[c(1L, 2L), c(1L, 2L)] <- TRUE
  out <- shift_mask(mat, vertical = -5L, horizontal = -5L)
  expect_equal(sum(out), 0L)
})

test_that("shift_mask validates inputs", {
  expect_error(shift_mask("not-a-mask"), "logical vector or logical matrix")
  m <- make_face_mask(c(16L, 16L), region = "full")
  expect_error(shift_mask(m), "img_dims")
  expect_error(
    shift_mask(m, vertical = 1L, img_dims = c(8L, 8L)),
    "length does not match"
  )
  expect_error(shift_mask(m, vertical = NA_real_, img_dims = c(16L, 16L)),
               "finite")
})

test_that("region_bounds_from_pixels round-trips against make_face_mask", {
  bounds <- region_bounds_from_pixels(
    row_min = 100, row_max = 130,
    col_min = 60,  col_max = 110,
    img_dims = c(256L, 256L)
  )
  expect_length(bounds, 4L)
  expect_true(all(bounds >= 0 & bounds <= 1))
  expect_lt(bounds[1L], bounds[2L])
  expect_lt(bounds[3L], bounds[4L])

  m <- make_face_mask(c(256L, 256L), region = "left_eye",
                      region_bounds = bounds)
  expect_true(any(m))
})

test_that("region_bounds_from_pixels rejects out-of-image or unordered input", {
  expect_error(
    region_bounds_from_pixels(0, 50, 60, 100, img_dims = 256L),
    "outside the image"
  )
  expect_error(
    region_bounds_from_pixels(50, 30, 60, 100, img_dims = 256L),
    "row_min <= row_max"
  )
  expect_error(
    region_bounds_from_pixels(50, 60, NA, 100, img_dims = 256L),
    "finite number"
  )
})
