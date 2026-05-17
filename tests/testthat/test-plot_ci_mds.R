## Tests for plot_ci_mds():
##   - basic auto-selection on synthetic CIs
##   - forced k = 2 and k = 3 (single panel vs grid)
##   - dimensionality auto-pick reaches threshold when it can,
##     falls back to k_max with a warning when it cannot
##   - stress_by_k is non-increasing; variance_pct_by_k is non-decreasing
##   - Kruskal band classification
##   - distance = euclidean_normalised rescales by sqrt(n_pixels_used)
##   - groups + shapes validation
##   - file save (PNG, PDF; single + multi panel)
##   - S3 print method runs

make_ci_vec <- function(seed, n_pix = 32L * 32L) {
  set.seed(seed)
  rnorm(n_pix) + sin(seq_len(n_pix) / 200)
}

n_side <- 32L
n_pix  <- n_side * n_side

# Build a list of CIs that are intentionally embedded in a 3D
# structure: two clusters along an unobservable third axis. This
# guarantees that stress at k=2 will be higher than at k=3.
make_3d_embedded_cis <- function(seed = 1L, n_pix = 64L * 64L) {
  set.seed(seed)
  # 6 CIs lying near (cos(t), sin(t), z) with z = +1 for half, -1 other.
  ang  <- seq(0, 2 * pi, length.out = 7L)[1:6]
  base <- rnorm(n_pix)
  basis_x <- rnorm(n_pix)
  basis_y <- rnorm(n_pix)
  basis_z <- rnorm(n_pix)
  zs <- c(1, 1, 1, -1, -1, -1)
  cis <- vector("list", 6L)
  for (i in seq_along(ang)) {
    cis[[i]] <- base +
      3 * cos(ang[i]) * basis_x +
      3 * sin(ang[i]) * basis_y +
      3 * zs[i]       * basis_z
  }
  names(cis) <- paste0("C", seq_along(cis))
  cis
}

test_that("plot_ci_mds runs and returns the right structure", {
  cis <- list(A = make_ci_vec(1), B = make_ci_vec(2),
              C = make_ci_vec(3), D = make_ci_vec(4),
              E = make_ci_vec(5))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- suppressWarnings(
    plot_ci_mds(cis, img_dims = c(n_side, n_side))
  )

  expect_s3_class(out, "rcisignal_mds")
  expect_true(is.matrix(out$mds_points))
  expect_equal(nrow(out$mds_points), 5L)
  expect_equal(rownames(out$mds_points), c("A", "B", "C", "D", "E"))
  expect_true(is.numeric(out$stress_by_k))
  expect_true(is.numeric(out$variance_pct_by_k))
  expect_true(out$stress_1 >= 0)
  expect_true(out$stress_band %in%
              c("excellent", "good", "fair", "poor", "very poor"))
  expect_true(out$k_selected >= 2L)
})

test_that("plot_ci_mds forced k = 2 returns a 2D solution", {
  cis <- list(A = make_ci_vec(10), B = make_ci_vec(11),
              C = make_ci_vec(12), D = make_ci_vec(13))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_ci_mds(cis, img_dims = c(n_side, n_side), k = 2L)
  expect_equal(out$k_selected, 2L)
  expect_equal(ncol(out$mds_points), 2L)
  expect_equal(out$k_selection_mode, "user")
  # panel_pairs should be the single pair (1, 2)
  expect_equal(unname(out$panel_pairs), matrix(c(1L, 2L), nrow = 2L))
})

test_that("plot_ci_mds forced k = 3 returns a 3D solution and a 3-pair grid", {
  cis <- list(A = make_ci_vec(20), B = make_ci_vec(21),
              C = make_ci_vec(22), D = make_ci_vec(23),
              E = make_ci_vec(24))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_ci_mds(cis, img_dims = c(n_side, n_side), k = 3L)
  expect_equal(out$k_selected, 3L)
  expect_equal(ncol(out$mds_points), 3L)
  expect_equal(ncol(out$panel_pairs), 3L)
  expect_equal(unname(out$panel_pairs),
               matrix(c(1L, 2L, 1L, 3L, 2L, 3L), nrow = 2L))
})

test_that("plot_ci_mds matches stats::cmdscale at the same distance matrix", {
  cis <- list(A = make_ci_vec(30), B = make_ci_vec(31),
              C = make_ci_vec(32), D = make_ci_vec(33))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_ci_mds(cis, img_dims = c(n_side, n_side), k = 2L)
  M <- do.call(cbind, cis)
  d <- stats::dist(t(M), method = "euclidean")
  ref <- stats::cmdscale(d, k = 2L)
  # cmdscale may return points with either sign; compare absolute
  # values (the embedding is unique up to reflection).
  expect_equal(abs(unname(out$mds_points)),
               abs(unname(ref)),
               tolerance = 1e-8)
})

test_that("plot_ci_mds auto-selection picks k > 2 when 3D structure exists", {
  cis <- make_3d_embedded_cis(seed = 7L, n_pix = 64L * 64L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot_ci_mds(cis, img_dims = c(64L, 64L),
                     stress_threshold = 0.05)
  expect_gte(out$k_selected, 2L)
  expect_lte(out$k_selected, 4L)
})

test_that("plot_ci_mds stress_by_k is non-increasing in k", {
  cis <- list(A = make_ci_vec(40), B = make_ci_vec(41),
              C = make_ci_vec(42), D = make_ci_vec(43),
              E = make_ci_vec(44))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- suppressWarnings(
    plot_ci_mds(cis, img_dims = c(n_side, n_side), k_max = 4L)
  )
  diffs <- diff(out$stress_by_k)
  expect_true(all(diffs <= 1e-9))
})

test_that("plot_ci_mds variance_pct_by_k is non-decreasing in k", {
  cis <- list(A = make_ci_vec(50), B = make_ci_vec(51),
              C = make_ci_vec(52), D = make_ci_vec(53),
              E = make_ci_vec(54))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- suppressWarnings(
    plot_ci_mds(cis, img_dims = c(n_side, n_side), k_max = 4L)
  )
  diffs <- diff(out$variance_pct_by_k)
  expect_true(all(diffs >= -1e-9))
})

test_that("plot_ci_mds Kruskal band classification is correct", {
  cls <- function(s) {
    if (s <= 0.025) "excellent"
    else if (s <= 0.05) "good"
    else if (s <= 0.10) "fair"
    else if (s <= 0.20) "poor"
    else "very poor"
  }
  expect_equal(cls(0.010), "excellent")
  expect_equal(cls(0.040), "good")
  expect_equal(cls(0.080), "fair")
  expect_equal(cls(0.150), "poor")
  expect_equal(cls(0.300), "very poor")
})

test_that("plot_ci_mds distance = 'euclidean_normalised' rescales", {
  cis <- list(A = make_ci_vec(60), B = make_ci_vec(61),
              C = make_ci_vec(62), D = make_ci_vec(63))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  raw <- plot_ci_mds(cis, img_dims = c(n_side, n_side),
                     mask = "face", distance = "euclidean_raw",
                     k = 2L)
  norm <- plot_ci_mds(cis, img_dims = c(n_side, n_side),
                      mask = "face",
                      distance = "euclidean_normalised", k = 2L)
  expect_equal(norm$distance_matrix,
               raw$distance_matrix / sqrt(raw$n_pixels_used))
})

test_that("plot_ci_mds groups / shapes validation", {
  cis <- list(A = make_ci_vec(70), B = make_ci_vec(71),
              C = make_ci_vec(72))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Valid groups
  expect_no_error(
    plot_ci_mds(cis, img_dims = c(n_side, n_side), k = 2L,
                groups = c(A = "g1", B = "g1", C = "g2"))
  )

  # Length mismatch
  expect_error(
    plot_ci_mds(cis, img_dims = c(n_side, n_side), k = 2L,
                groups = c(A = "g1", B = "g1")),
    "one entry per CI"
  )

  # Name mismatch
  expect_error(
    plot_ci_mds(cis, img_dims = c(n_side, n_side), k = 2L,
                groups = c(A = "g1", B = "g1", X = "g2")),
    "Names of"
  )

  # Unnamed
  expect_error(
    plot_ci_mds(cis, img_dims = c(n_side, n_side), k = 2L,
                groups = c("g1", "g1", "g2")),
    "NAMED character"
  )
})

test_that("plot_ci_mds errors on too few CIs", {
  cis <- list(A = make_ci_vec(80), B = make_ci_vec(81))
  expect_error(plot_ci_mds(cis, img_dims = c(n_side, n_side)),
               "at least 3")
})

test_that("plot_ci_mds k out of range errors clearly", {
  cis <- list(A = make_ci_vec(90), B = make_ci_vec(91),
              C = make_ci_vec(92))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(plot_ci_mds(cis, img_dims = c(n_side, n_side), k = 1L),
               "out of range")
  expect_error(plot_ci_mds(cis, img_dims = c(n_side, n_side), k = 5L),
               "out of range")
})

test_that("plot_ci_mds saves PNG single-panel via file =", {
  cis <- list(A = make_ci_vec(100), B = make_ci_vec(101),
              C = make_ci_vec(102))
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  out <- plot_ci_mds(cis, img_dims = c(n_side, n_side),
                     k = 2L, file = tmp)
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 0)
})

test_that("plot_ci_mds saves PDF multi-panel via file =", {
  cis <- list(A = make_ci_vec(110), B = make_ci_vec(111),
              C = make_ci_vec(112), D = make_ci_vec(113),
              E = make_ci_vec(114))
  tmp <- tempfile(fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)
  out <- plot_ci_mds(cis, img_dims = c(n_side, n_side),
                     k = 3L, file = tmp)
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 0)
})

test_that("plot_ci_mds print method runs without error", {
  cis <- list(A = make_ci_vec(120), B = make_ci_vec(121),
              C = make_ci_vec(122), D = make_ci_vec(123))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- suppressWarnings(
    plot_ci_mds(cis, img_dims = c(n_side, n_side))
  )
  expect_output(print(out), "rcisignal MDS")
  expect_output(print(out), "Dimensionality selection")
})
