test_that("group_ci returns one column per group with correct means", {
  sm <- make_sig(n_pix = 256L, n_p = 8L, seed = 1L)
  g  <- rep(c("A", "B"), each = 4L)
  gcis <- group_ci(sm, by = g)

  expect_s3_class(gcis, "rcisignal_group_ci")
  expect_true(is.matrix(gcis))
  expect_equal(dim(gcis), c(256L, 2L))
  expect_equal(colnames(gcis), c("A", "B"))

  expect_equal(as.numeric(gcis[, "A"]),
               as.numeric(rowMeans(sm[, 1:4])))
  expect_equal(as.numeric(gcis[, "B"]),
               as.numeric(rowMeans(sm[, 5:8])))
  expect_identical(attr(gcis, "n"),
                   c(A = 4L, B = 4L))
})

test_that("group_ci factorial grouping joins level names with _", {
  sm <- make_sig(n_pix = 64L, n_p = 8L, seed = 2L)
  country <- rep(c("US", "PT"), each = 4L)
  q       <- rep(c("q1", "q2"), times = 4L)
  gcis <- group_ci(sm, by = list(country = country, q = q))
  expect_setequal(colnames(gcis),
                  c("US_q1", "US_q2", "PT_q1", "PT_q2"))
  expect_equal(sum(attr(gcis, "n")), 8L)
})

test_that("group_ci drop = FALSE keeps empty cells as NA columns", {
  sm <- make_sig(n_pix = 64L, n_p = 4L, seed = 3L)
  g  <- factor(c("A", "A", "B", "B"), levels = c("A", "B", "C"))
  gcis_drop <- group_ci(sm, by = g, drop = TRUE)
  expect_equal(colnames(gcis_drop), c("A", "B"))

  gcis_keep <- group_ci(sm, by = g, drop = FALSE)
  expect_equal(colnames(gcis_keep), c("A", "B", "C"))
  expect_true(all(is.na(gcis_keep[, "C"])))
  expect_identical(attr(gcis_keep, "n")[["C"]], 0L)
})

test_that("group_ci aborts on length mismatch with clear message", {
  sm <- make_sig(n_pix = 64L, n_p = 6L, seed = 4L)
  expect_error(group_ci(sm, by = c("A", "B")),
               regexp = "length")
})

test_that("group_ci warns and drops producers with NA group", {
  sm <- make_sig(n_pix = 64L, n_p = 6L, seed = 5L)
  g  <- c("A", "A", "A", "B", "B", NA)
  expect_warning(
    gcis <- group_ci(sm, by = g),
    regexp = "NA"
  )
  expect_setequal(colnames(gcis), c("A", "B"))
  expect_equal(attr(gcis, "n")[["A"]], 3L)
  expect_equal(attr(gcis, "n")[["B"]], 2L)
})

test_that("group_ci preserves img_dims attribute", {
  sm <- make_sig(n_pix = 256L, n_p = 4L, seed = 6L)
  gcis <- group_ci(sm, by = c("A", "A", "B", "B"))
  expect_equal(attr(gcis, "img_dims"),
               attr(sm, "img_dims"))
})

test_that("group_ci aborts when handed its own output", {
  sm <- make_sig(n_pix = 64L, n_p = 4L, seed = 7L)
  gcis <- group_ci(sm, by = c("A", "A", "B", "B"))
  expect_error(group_ci(gcis, by = c("X", "Y")),
               regexp = "stage-2")
})

test_that("print.rcisignal_group_ci surfaces group sizes", {
  sm <- make_sig(n_pix = 64L, n_p = 6L, seed = 8L)
  gcis <- group_ci(sm, by = c("A", "A", "B", "B", "C", "C"))
  out <- utils::capture.output(print(gcis))
  expect_true(any(grepl("3 groups", out)))
  expect_true(any(grepl("\\bA\\b\\s+n = 2", out)))
})

test_that("as.list.rcisignal_group_ci returns named per-group vectors", {
  sm <- make_sig(n_pix = 64L, n_p = 4L, seed = 9L)
  gcis <- group_ci(sm, by = c("A", "A", "B", "B"))
  L <- as.list(gcis)
  expect_named(L, c("A", "B"))
  expect_equal(length(L[[1L]]), nrow(sm))
})

test_that("group_ci output feeds plot_ci_distance_matrix", {
  sm <- make_sig(n_pix = 256L, n_p = 8L, seed = 10L)
  gcis <- group_ci(sm, by = rep(c("A", "B", "C", "D"), each = 2L))
  out <- plot_ci_distance_matrix(gcis,
                                 img_dims    = c(16L, 16L),
                                 show_values = FALSE,
                                 file        = tempfile(fileext = ".pdf"))
  expect_type(out, "list")
  expect_true(is.matrix(out$distance_matrix))
  expect_equal(dim(out$distance_matrix), c(4L, 4L))
})
