## save_ci_images() round-trip + filename convention.

skip_if_not_installed("png")
skip_if_not_installed("jpeg")

make_ci_set <- function(n_pix = 1024L, n_p = 4L,
                        level = c("individual", "group"),
                        seed = 1L) {
  level <- match.arg(level)
  sm <- make_sig(n_pix = n_pix, n_p = n_p, seed = seed)
  colnames(sm) <- if (level == "individual") {
    sprintf("P%03d", seq_len(n_p))
  } else {
    sprintf("cond_%d", seq_len(n_p))
  }
  attr(sm, "ci_level") <- level
  sm
}

make_base <- function(img_dims) {
  matrix(0.5, img_dims[1L], img_dims[2L])
}

test_that("save_ci_images writes one PNG per column with ind_ci_ prefix", {
  sm <- make_ci_set(n_pix = 1024L, n_p = 3L, level = "individual")
  out <- tempfile("ci_test_")
  base <- make_base(attr(sm, "img_dims"))

  paths <- save_ci_images(sm, base, out, quiet = TRUE)
  expect_length(paths, 3L)
  expect_true(all(file.exists(paths)))
  expect_match(basename(paths), "^ind_ci_P\\d{3}\\.png$", all = TRUE)

  img <- png::readPNG(paths[1L])
  # Default palette is "grayscale": writePNG of a 2D matrix returns a
  # 2D matrix on read (no third color-channel dim).
  expect_length(dim(img), 2L)
  expect_equal(dim(img), attr(sm, "img_dims"))
})

test_that("grayscale default matches rcicr's applyScaling + combine pipeline", {
  # Replicate rcicr::generateCI() PNG output exactly for one CI:
  # scaled = applyScaling(base, ci, "independent", 0.1)
  # combined = (scaled + base) / 2
  # png::writePNG(combined)
  sm <- make_ci_set(n_pix = 1024L, n_p = 1L, level = "individual")
  img_dims <- attr(sm, "img_dims")
  set.seed(42L)
  base <- matrix(runif(prod(img_dims), 0.3, 0.7),
                 img_dims[1L], img_dims[2L])

  out <- tempfile("rcicr_match_"); dir.create(out)
  paths <- save_ci_images(sm, base, out, quiet = TRUE)
  written <- png::readPNG(paths[1L])

  # Reference: rcicr's exact pipeline (independent scaling).
  ci <- as.vector(sm[, 1L])
  k <- max(abs(range(ci)))
  scaled <- (ci + k) / (2 * k)
  combined <- (scaled + as.vector(base)) / 2
  combined[combined < 0] <- 0
  combined[combined > 1] <- 1
  reference <- matrix(combined, img_dims[1L], img_dims[2L])

  # 8-bit PNG quantization tolerance: 1/255.
  expect_equal(written, reference, tolerance = 1.5 / 255)
})

test_that("grayscale + scaling = 'constant' matches rcicr formula", {
  sm <- make_ci_set(n_pix = 1024L, n_p = 1L, level = "individual")
  img_dims <- attr(sm, "img_dims")
  base <- matrix(0.5, img_dims[1L], img_dims[2L])

  out <- tempfile("rcicr_const_"); dir.create(out)
  paths <- suppressWarnings(
    save_ci_images(sm, base, out, scaling = "constant",
                   scaling_constant = 0.5, quiet = TRUE)
  )
  written <- png::readPNG(paths[1L])

  ci <- as.vector(sm[, 1L])
  scaled <- (ci + 0.5) / (2 * 0.5)
  combined <- (scaled + as.vector(base)) / 2
  combined[combined < 0] <- 0
  combined[combined > 1] <- 1
  reference <- matrix(combined, img_dims[1L], img_dims[2L])

  expect_equal(written, reference, tolerance = 1.5 / 255)
})

test_that("grayscale + scaling = 'constant' warns when out of [0,1]", {
  sm <- make_ci_set(n_pix = 1024L, n_p = 1L, level = "individual",
                    seed = 7L)
  # Inflate so (ci + 0.01) / 0.02 blows past 1.
  sm[, 1L] <- sm[, 1L] * 10
  img_dims <- attr(sm, "img_dims")
  base <- matrix(0.5, img_dims[1L], img_dims[2L])
  out <- tempfile("rcicr_warn_"); dir.create(out)
  expect_warning(
    save_ci_images(sm, base, out, scaling = "constant",
                   scaling_constant = 0.01, quiet = TRUE),
    regexp = "clipping|outside"
  )
})

test_that("palette = 'diverging' still produces an RGB image", {
  sm <- make_ci_set(n_pix = 1024L, n_p = 1L, level = "individual")
  base <- make_base(attr(sm, "img_dims"))
  out <- tempfile("ci_div_")
  paths <- save_ci_images(sm, base, out, palette = "diverging",
                          quiet = TRUE)
  img <- png::readPNG(paths[1L])
  expect_length(dim(img), 3L)
  expect_equal(dim(img)[1:2], attr(sm, "img_dims"))
})

test_that("save_ci_images uses group_ci_ prefix on a group-level matrix", {
  sm <- make_ci_set(n_pix = 1024L, n_p = 2L, level = "group")
  out <- tempfile("ci_test_grp_")
  base <- make_base(attr(sm, "img_dims"))

  paths <- save_ci_images(sm, base, out, quiet = TRUE)
  expect_match(basename(paths), "^group_ci_cond_\\d\\.png$", all = TRUE)
})

test_that("custom prefix overrides level detection", {
  sm <- make_ci_set(level = "group")
  out <- tempfile("ci_test_pref_")
  base <- make_base(attr(sm, "img_dims"))
  paths <- save_ci_images(sm, base, out, prefix = "trust_",
                          quiet = TRUE)
  expect_match(basename(paths), "^trust_", all = TRUE)
})

test_that("jpeg format works and accepts quality", {
  sm <- make_ci_set(n_p = 2L)
  out <- tempfile("ci_test_jpg_")
  base <- make_base(attr(sm, "img_dims"))
  paths <- save_ci_images(sm, base, out, format = "jpeg",
                          quality = 80, quiet = TRUE)
  expect_true(all(file.exists(paths)))
  expect_match(basename(paths), "\\.jpeg$", all = TRUE)
})

test_that("save_ci_images aborts on duplicate or missing colnames", {
  sm <- make_ci_set(n_p = 3L)
  base <- make_base(attr(sm, "img_dims"))
  colnames(sm) <- c("A", "A", "B")
  expect_error(save_ci_images(sm, base, tempfile("dup_")),
               regexp = "duplicat")

  sm2 <- make_ci_set(n_p = 3L)
  colnames(sm2) <- NULL
  expect_error(save_ci_images(sm2, base, tempfile("no_names_")),
               regexp = "column names")
})

test_that("save_ci_images refuses overwrite by default", {
  sm <- make_ci_set(n_p = 2L)
  out <- tempfile("ci_test_ovw_")
  base <- make_base(attr(sm, "img_dims"))

  save_ci_images(sm, base, out, quiet = TRUE)
  expect_error(save_ci_images(sm, base, out, quiet = TRUE),
               regexp = "exist|overwrite")

  # overwrite = TRUE silences it
  paths <- save_ci_images(sm, base, out, overwrite = TRUE,
                          quiet = TRUE)
  expect_true(all(file.exists(paths)))
})

test_that("save_ci_images creates the dir if missing", {
  sm <- make_ci_set(n_p = 1L)
  base <- make_base(attr(sm, "img_dims"))
  nested <- file.path(tempfile("mk_"), "deeper", "dir")
  paths <- save_ci_images(sm, base, nested, quiet = TRUE)
  expect_true(dir.exists(nested))
  expect_true(file.exists(paths))
})

test_that("save_ci_images works with the fire palette", {
  sm <- make_ci_set(n_p = 1L)
  base <- make_base(attr(sm, "img_dims"))
  out <- tempfile("ci_fire_")
  paths <- save_ci_images(sm, base, out, palette = "fire",
                          quiet = TRUE)
  img <- png::readPNG(paths)
  # YlOrRd ramp; red channel should dominate green at high values
  expect_true(mean(img[, , 1]) >= mean(img[, , 2]))
})
