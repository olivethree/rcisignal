test_that("Brief-RC mask equals Schmitz genMask by hand", {
  skip_if_not_installed("png")
  set.seed(1)
  n_pix <- 16L * 16L
  n_pool <- 100L
  noise_matrix <- matrix(rnorm(n_pix * n_pool), n_pix, n_pool)

  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 16L, 16L), base_path)

  responses <- data.frame(
    participant_id = rep(c("p1", "p2"), each = 30L),
    trial          = rep(1:30, 2L),
    stimulus       = sample.int(n_pool, 60L, replace = TRUE),
    response       = sample(c(-1L, 1L), 60L, replace = TRUE)
  )

  res <- ci_from_responses_briefrc(
    responses       = responses,
    noise_matrix    = noise_matrix,
    base_image      = base_path
  )
  expect_equal(dim(res$signal_matrix), c(n_pix, 2L))

  # hand computation for p1
  sub <- responses[responses$participant_id == "p1", ]
  X <- data.table::data.table(response = sub$response, stim = sub$stimulus)
  X <- X[, list(response = mean(response)), by = "stim"]
  hand_mask <- as.numeric(noise_matrix[, X$stim] %*% X$response) / nrow(X)
  expect_equal(res$signal_matrix[, "p1"], hand_mask, tolerance = 1e-10)
})

test_that("duplicate-stimulus responses average (genMask rule)", {
  skip_if_not_installed("png")
  n_pix <- 64L
  noise_matrix <- matrix(c(rep(1, n_pix), rep(2, n_pix), rep(3, n_pix)),
                         nrow = n_pix, ncol = 3L)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 8L, 8L), base_path)

  # p1 chose stim=1 twice with opposite responses (should cancel)
  # and stim=2 once with response=1 (unique).
  # Expected: mean(response) by stim -> stim=1 weight 0, stim=2 weight 1
  # mask = noise[,1]*0 + noise[,2]*1, divided by 2 unique stims
  responses <- data.frame(
    participant_id = "p1",
    trial          = 1:3,
    stimulus       = c(1L, 1L, 2L),
    response       = c(1L, -1L, 1L)
  )
  res <- ci_from_responses_briefrc(
    responses       = responses,
    noise_matrix    = noise_matrix,
    base_image      = base_path
  )
  expected <- (noise_matrix[, 1] * 0 + noise_matrix[, 2] * 1) / 2
  expect_equal(as.numeric(res$signal_matrix[, "p1"]), expected,
               tolerance = 1e-10)
})

test_that("Brief-RC rejects bad response coding", {
  skip_if_not_installed("png")
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 4L, 4L), base_path)
  noise_matrix <- matrix(0, nrow = 16L, ncol = 5L)
  responses <- data.frame(
    participant_id = "p1", trial = 1:3,
    stimulus = c(1L, 2L, 3L),
    response = c(0L, 1L, -1L)   # 0 is invalid
  )
  expect_error(
    ci_from_responses_briefrc(responses, noise_matrix = noise_matrix,
                              base_image      = base_path),
    "must contain only values"
  )
})

test_that("Brief-RC `{0,1}` miscoding gets a recoding hint in error", {
  skip_if_not_installed("png")
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 4L, 4L), base_path)
  noise_matrix <- matrix(0, nrow = 16L, ncol = 5L)
  responses <- data.frame(
    participant_id = rep("p1", 4L), trial = 1:4,
    stimulus = c(1L, 2L, 3L, 4L),
    response = c(0L, 1L, 0L, 1L)   # the classic miscoding
  )
  err <- tryCatch(
    ci_from_responses_briefrc(responses, noise_matrix = noise_matrix,
                              base_image      = base_path),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "Did you mean")
  expect_match(err, "Recode in one line")
})

test_that("briefrc20 produces same signal matrix as briefrc12 on identical data", {
  # Mathematically the formula is symmetric in the per-trial split
  # (6/6 for briefrc12, 10/10 for briefrc20). The `method`
  # argument is metadata: it should not branch the computation.
  skip_if_not_installed("png")
  set.seed(11)
  n_pix <- 16L * 16L
  noise_matrix <- matrix(rnorm(n_pix * 50L), n_pix, 50L)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 16L, 16L), base_path)
  responses <- data.frame(
    participant_id = rep(c("p1", "p2"), each = 30L),
    trial          = rep(1:30, 2L),
    stimulus       = sample.int(50L, 60L, replace = TRUE),
    response       = sample(c(-1L, 1L), 60L, replace = TRUE)
  )

  res12 <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix,
    base_image      = base_path, method = "briefrc12"
  )
  res20 <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix,
    base_image      = base_path, method = "briefrc20"
  )
  expect_equal(res20$signal_matrix, res12$signal_matrix)
  expect_identical(res20$method, "briefrc20")
  expect_identical(res12$method, "briefrc12")
})

test_that("briefrc20 round-trips: hand-computed mask matches genMask", {
  # Positive round-trip: feed a 20-image-per-trial-shaped dataset
  # and verify per-producer mask equals the Schmitz formula
  # applied by hand.
  skip_if_not_installed("png")
  set.seed(12)
  n_pix <- 16L * 16L
  n_pool <- 200L
  noise_matrix <- matrix(rnorm(n_pix * n_pool), n_pix, n_pool)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 16L, 16L), base_path)

  # Brief-RC 20: 50 trials per producer; on each trial the
  # producer picks one stimulus from the 20-image set and labels
  # it +1 (original) or -1 (inverted). Recorded as one row per
  # trial with the chosen pool id.
  responses <- data.frame(
    participant_id = rep(c("p1", "p2"), each = 50L),
    trial          = rep(1:50, 2L),
    stimulus       = sample.int(n_pool, 100L, replace = TRUE),
    response       = sample(c(-1L, 1L), 100L, replace = TRUE)
  )

  res <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix,
    base_image      = base_path, method = "briefrc20"
  )
  expect_equal(dim(res$signal_matrix), c(n_pix, 2L))

  sub <- responses[responses$participant_id == "p1", ]
  X <- data.table::data.table(response = sub$response, stim = sub$stimulus)
  X <- X[, list(response = mean(response)), by = "stim"]
  hand_mask <- as.numeric(noise_matrix[, X$stim] %*% X$response) / nrow(X)
  expect_equal(res$signal_matrix[, "p1"], hand_mask, tolerance = 1e-10)
})

test_that("default scaling = 'none' returns no rendered_ci field", {
  skip_if_not_installed("png")
  set.seed(2)
  n_pix <- 8L * 8L
  noise_matrix <- matrix(rnorm(n_pix * 20L), n_pix, 20L)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 8L, 8L), base_path)
  responses <- data.frame(
    participant_id = rep("p1", 10L),
    trial          = 1:10,
    stimulus       = sample.int(20L, 10L),
    response       = sample(c(-1L, 1L), 10L, replace = TRUE)
  )
  res <- ci_from_responses_briefrc(
    responses       = responses,
    noise_matrix    = noise_matrix,
    base_image      = base_path
  )
  expect_null(res$rendered_ci)
  expect_equal(res$scaling, "none")
})

test_that("scaling = 'matched' adds rendered_ci, signal_matrix unchanged", {
  skip_if_not_installed("png")
  set.seed(3)
  n_pix <- 8L * 8L
  noise_matrix <- matrix(rnorm(n_pix * 20L), n_pix, 20L)
  # base must have within-image variance for "matched" to do anything
  base_img <- matrix(seq(0, 1, length.out = n_pix), 8L, 8L)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(base_img, base_path)
  responses <- data.frame(
    participant_id = rep(c("p1", "p2"), each = 8L),
    trial          = c(1:8, 1:8),
    stimulus       = sample.int(20L, 16L, replace = TRUE),
    response       = sample(c(-1L, 1L), 16L, replace = TRUE)
  )
  res_none <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix,
    base_image      = base_path, scaling = "none"
  )
  res_match <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix,
    base_image      = base_path, scaling = "matched"
  )
  expect_equal(res_match$signal_matrix, res_none$signal_matrix)
  expect_false(is.null(res_match$rendered_ci))
  expect_equal(dim(res_match$rendered_ci), dim(res_match$signal_matrix))
  # matched stretches mask range to base range, then adds base.
  # Per-column rendered range should be at least the base range
  # (since rendered = base + matched(mask), and matched(mask) has
  # the same range as base).
  base_rng <- diff(range(base_img))
  rendered_rngs <- vapply(seq_len(ncol(res_match$rendered_ci)),
                          function(j) diff(range(res_match$rendered_ci[, j])),
                          numeric(1L))
  expect_true(all(rendered_rngs >= base_rng - 1e-3))
})

test_that("scaling = 'constant' applies the requested multiplier", {
  skip_if_not_installed("png")
  set.seed(4)
  n_pix <- 8L * 8L
  noise_matrix <- matrix(rnorm(n_pix * 20L), n_pix, 20L)
  # Use a base value that round-trips PNG losslessly: png stores
  # 8-bit integers so 0.5 becomes 127/255 != 0.5.
  base_img <- matrix(127 / 255, 8L, 8L)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(base_img, base_path)
  responses <- data.frame(
    participant_id = rep("p1", 10L),
    trial          = 1:10,
    stimulus       = sample.int(20L, 10L),
    response       = sample(c(-1L, 1L), 10L, replace = TRUE)
  )
  res <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix,
    base_image       = base_path,
    scaling          = "constant",
    scaling_constant = 3
  )
  # rendered = base_read_back + 3 * mask
  base_read <- as.vector(png::readPNG(base_path))
  expect_equal(res$rendered_ci[, 1L],
               base_read + 3 * res$signal_matrix[, 1L],
               tolerance = 1e-6)
})

test_that("scaling = 'constant' without scaling_constant aborts", {
  skip_if_not_installed("png")
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 4L, 4L), base_path)
  responses <- data.frame(
    participant_id = "p1", trial = 1L, stimulus = 1L, response = 1L
  )
  expect_error(
    ci_from_responses_briefrc(
      responses, noise_matrix = matrix(0, 16L, 2L),
      base_image      = base_path,
      scaling = "constant"
    ),
    "scaling_constant"
  )
})

test_that("base_image accepts a numeric matrix", {
  set.seed(21)
  n_pix <- 16L * 16L
  noise_matrix <- matrix(rnorm(n_pix * 50L), n_pix, 50L)
  base_mat <- matrix(seq(0, 1, length.out = n_pix), 16L, 16L)
  responses <- data.frame(
    participant_id = rep(c("p1", "p2"), each = 30L),
    trial          = rep(1:30, 2L),
    stimulus       = sample.int(50L, 60L, replace = TRUE),
    response       = sample(c(-1L, 1L), 60L, replace = TRUE)
  )
  res <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix,
    base_image = base_mat, scaling = "matched"
  )
  expect_equal(dim(res$signal_matrix), c(n_pix, 2L))
  expect_equal(res$img_dims, c(16L, 16L))
  expect_false(is.null(res$rendered_ci))
})

test_that("base_image matrix and equivalent PNG path agree", {
  skip_if_not_installed("png")
  set.seed(22)
  n_pix <- 16L * 16L
  noise_matrix <- matrix(rnorm(n_pix * 50L), n_pix, 50L)
  base_mat <- matrix(127 / 255, 16L, 16L)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(base_mat, base_path)
  responses <- data.frame(
    participant_id = rep("p1", 20L),
    trial          = 1:20,
    stimulus       = sample.int(50L, 20L),
    response       = sample(c(-1L, 1L), 20L, replace = TRUE)
  )
  res_mat  <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix, base_image = base_mat
  )
  res_path <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix, base_image = base_path
  )
  expect_equal(res_mat$signal_matrix, res_path$signal_matrix,
               tolerance = 1e-10)
})

test_that("base_image is optional when scaling = 'none'", {
  set.seed(23)
  n_pix <- 16L * 16L
  noise_matrix <- matrix(rnorm(n_pix * 30L), n_pix, 30L)
  responses <- data.frame(
    participant_id = rep("p1", 10L),
    trial          = 1:10,
    stimulus       = sample.int(30L, 10L),
    response       = sample(c(-1L, 1L), 10L, replace = TRUE)
  )
  res <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix
  )
  expect_equal(res$img_dims, c(16L, 16L))
  expect_null(res$rendered_ci)
})

test_that("missing base_image with scaling != 'none' aborts", {
  noise_matrix <- matrix(0, nrow = 16L, ncol = 5L)
  responses <- data.frame(
    participant_id = "p1",
    trial          = 1:2,
    stimulus       = c(1L, 2L),
    response       = c(1L, -1L)
  )
  expect_error(
    ci_from_responses_briefrc(
      responses, noise_matrix = noise_matrix, scaling = "matched"
    ),
    "base_image"
  )
})

test_that("base_image rejects out-of-range matrix values", {
  noise_matrix <- matrix(0, nrow = 16L, ncol = 5L)
  responses <- data.frame(
    participant_id = "p1",
    trial          = 1:2,
    stimulus       = c(1L, 2L),
    response       = c(1L, -1L)
  )
  bad_base <- matrix(2, 4L, 4L)
  expect_error(
    ci_from_responses_briefrc(
      responses, noise_matrix = noise_matrix, base_image = bad_base
    ),
    "\\[0, 1\\]"
  )
})

test_that("base_image_path is a deprecated alias for base_image", {
  skip_if_not_installed("png")
  set.seed(24)
  n_pix <- 16L * 16L
  noise_matrix <- matrix(rnorm(n_pix * 30L), n_pix, 30L)
  base_path <- tempfile(fileext = ".png")
  png::writePNG(matrix(0.5, 16L, 16L), base_path)
  responses <- data.frame(
    participant_id = rep("p1", 10L),
    trial          = 1:10,
    stimulus       = sample.int(30L, 10L),
    response       = sample(c(-1L, 1L), 10L, replace = TRUE)
  )
  res_old <- suppressWarnings(
    ci_from_responses_briefrc(
      responses, noise_matrix = noise_matrix,
      base_image_path = base_path
    )
  )
  expect_warning(
    ci_from_responses_briefrc(
      responses, noise_matrix = noise_matrix,
      base_image_path = base_path
    ),
    "deprecated"
  )
  res_new <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix, base_image = base_path
  )
  expect_equal(res_old$signal_matrix, res_new$signal_matrix)
})

test_that("group_by= returns both signal_matrix and group_ci", {
  set.seed(1)
  n_pix <- 16L * 16L
  noise_matrix <- matrix(rnorm(n_pix * 100L), n_pix, 100L)
  responses <- data.frame(
    participant_id = rep(c("p1", "p2", "p3", "p4"), each = 20L),
    condition      = rep(c("A", "A", "B", "B"), each = 20L),
    trial          = rep(seq_len(20L), 4L),
    stimulus       = sample.int(100L, 80L, replace = TRUE),
    response       = sample(c(-1L, 1L), 80L, replace = TRUE)
  )
  res <- ci_from_responses_briefrc(
    responses, noise_matrix = noise_matrix, group_by = "condition"
  )
  expect_true(!is.null(res$signal_matrix))
  expect_true(!is.null(res$group_ci))
  expect_equal(ncol(res$signal_matrix), 4L)
  expect_equal(ncol(res$group_ci), 2L)
  expect_setequal(colnames(res$group_ci), c("A", "B"))
  expect_identical(attr(res$signal_matrix, "ci_level"), "individual")
  expect_identical(attr(res$group_ci, "ci_level"), "group")
})
