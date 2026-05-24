## group_ci() v0.2.0 API: data frame + column-name idiom.

make_sig_with_pids <- function(n_pix = 64L, n_p = 8L, seed = 1L) {
  sm <- make_sig(n_pix = n_pix, n_p = n_p, seed = seed)
  colnames(sm) <- sprintf("p%02d", seq_len(n_p))
  sm
}

make_responses_two_conditions <- function(n_p = 8L, n_trials = 5L,
                                          seed = 1L) {
  pids <- sprintf("p%02d", seq_len(n_p))
  cond <- rep(c("A", "B"), each = n_p / 2L)
  set.seed(seed)
  data.frame(
    participant_id = rep(pids, each = n_trials),
    condition      = rep(cond, each = n_trials),
    response       = sample(c(-1L, 1L),
                            n_p * n_trials, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

test_that("group_ci returns one column per group with correct means", {
  sm <- make_sig_with_pids(n_pix = 256L, n_p = 8L, seed = 1L)
  responses <- make_responses_two_conditions(n_p = 8L)
  gcis <- group_ci(sm, responses, by = "condition")

  expect_s3_class(gcis, "rcisignal_group_ci")
  expect_true(is.matrix(gcis))
  expect_equal(dim(gcis), c(256L, 2L))
  expect_setequal(colnames(gcis), c("A", "B"))

  expect_equal(as.numeric(gcis[, "A"]),
               as.numeric(rowMeans(sm[, 1:4])))
  expect_equal(as.numeric(gcis[, "B"]),
               as.numeric(rowMeans(sm[, 5:8])))
  expect_identical(attr(gcis, "n"),
                   c(A = 4L, B = 4L))
})

test_that("group_ci factorial grouping joins levels with _", {
  sm <- make_sig_with_pids(n_pix = 64L, n_p = 8L, seed = 2L)
  # 8 producers, 2 trials each. Each producer is fully in one
  # country x q cell so the "consistent by within producer" check
  # is satisfied.
  responses <- data.frame(
    participant_id = rep(sprintf("p%02d", 1:8), each = 2L),
    country        = rep(c("US", "US", "US", "US", "PT", "PT", "PT", "PT"),
                         each = 2L),
    q              = rep(c("q1", "q1", "q2", "q2", "q1", "q1", "q2", "q2"),
                         each = 2L),
    stringsAsFactors = FALSE
  )
  gcis <- group_ci(sm, responses, by = c("country", "q"))
  expect_setequal(colnames(gcis),
                  c("US_q1", "US_q2", "PT_q1", "PT_q2"))
  expect_equal(sum(attr(gcis, "n")), 8L)
  expect_match(attr(gcis, "by_name"), "country x q", fixed = TRUE)
})

test_that("group_ci drop = FALSE keeps empty cells as NA columns", {
  sm <- make_sig_with_pids(n_pix = 64L, n_p = 4L, seed = 3L)
  responses <- data.frame(
    participant_id = rep(sprintf("p%02d", 1:4), each = 3L),
    condition      = rep(c("A", "A", "B", "B"), each = 3L),
    stringsAsFactors = FALSE
  )
  responses$condition <- factor(responses$condition,
                                levels = c("A", "B", "C"))

  gcis_drop <- group_ci(sm, responses, by = "condition", drop = TRUE)
  expect_setequal(colnames(gcis_drop), c("A", "B"))

  gcis_keep <- group_ci(sm, responses, by = "condition", drop = FALSE)
  expect_setequal(colnames(gcis_keep), c("A", "B", "C"))
  expect_true(all(is.na(gcis_keep[, "C"])))
  expect_identical(attr(gcis_keep, "n")[["C"]], 0L)
})

test_that("group_ci aborts when a producer is missing from responses", {
  sm <- make_sig_with_pids(n_pix = 64L, n_p = 4L, seed = 4L)
  responses <- data.frame(
    participant_id = rep(c("p01", "p02", "p03"), each = 3L),
    condition      = rep(c("A", "A", "B"), each = 3L),
    stringsAsFactors = FALSE
  )
  expect_error(group_ci(sm, responses, by = "condition"),
               regexp = "p04|not found")
})

test_that("group_ci aborts on inconsistent by within a producer", {
  sm <- make_sig_with_pids(n_pix = 64L, n_p = 4L, seed = 5L)
  responses <- data.frame(
    participant_id = c("p01", "p01", "p02", "p02",
                       "p03", "p03", "p04", "p04"),
    condition      = c("A", "B",                # p01 spans two conds
                       "A", "A",
                       "B", "B",
                       "B", "B"),
    stringsAsFactors = FALSE
  )
  expect_error(group_ci(sm, responses, by = "condition"),
               regexp = "inconsistent")
})

test_that("group_ci aborts on bad by", {
  sm <- make_sig_with_pids(n_pix = 64L, n_p = 4L, seed = 6L)
  responses <- data.frame(
    participant_id = sprintf("p%02d", 1:4),
    condition      = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  expect_error(group_ci(sm, responses, by = "no_such_col"),
               regexp = "no_such_col|not in")
  expect_error(group_ci(sm, responses, by = character(0L)),
               regexp = "character vector")
})

test_that("group_ci aborts when signal_matrix has no colnames", {
  sm <- make_sig(n_pix = 64L, n_p = 4L, seed = 7L)  # no colnames set
  responses <- data.frame(
    participant_id = sprintf("p%02d", 1:4),
    condition      = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  expect_error(group_ci(sm, responses, by = "condition"),
               regexp = "column names")
})

test_that("group_ci honors a non-default col_participant", {
  sm <- make_sig_with_pids(n_pix = 64L, n_p = 4L, seed = 8L)
  responses <- data.frame(
    subject_id = sprintf("p%02d", 1:4),
    condition  = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  gcis <- group_ci(sm, responses, by = "condition",
                   col_participant = "subject_id")
  expect_setequal(colnames(gcis), c("A", "B"))
})

test_that("group_ci preserves img_dims attribute", {
  sm <- make_sig_with_pids(n_pix = 256L, n_p = 4L, seed = 9L)
  responses <- data.frame(
    participant_id = sprintf("p%02d", 1:4),
    condition      = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  gcis <- group_ci(sm, responses, by = "condition")
  expect_equal(attr(gcis, "img_dims"),
               attr(sm, "img_dims"))
})

test_that("group_ci aborts when handed its own output", {
  sm <- make_sig_with_pids(n_pix = 64L, n_p = 4L, seed = 10L)
  responses <- data.frame(
    participant_id = sprintf("p%02d", 1:4),
    condition      = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  gcis <- group_ci(sm, responses, by = "condition")
  expect_error(group_ci(gcis, responses, by = "condition"),
               regexp = "stage-2")
})

test_that("print.rcisignal_group_ci surfaces group sizes", {
  sm <- make_sig_with_pids(n_pix = 64L, n_p = 6L, seed = 11L)
  responses <- data.frame(
    participant_id = rep(sprintf("p%02d", 1:6), each = 4L),
    condition      = rep(c("A", "A", "B", "B", "C", "C"), each = 4L),
    stringsAsFactors = FALSE
  )
  gcis <- group_ci(sm, responses, by = "condition")
  out <- utils::capture.output(print(gcis))
  expect_true(any(grepl("3 groups", out)))
  expect_true(any(grepl("\\bA\\b\\s+n = 2", out)))
})

test_that("as.list.rcisignal_group_ci returns named per-group vectors", {
  sm <- make_sig_with_pids(n_pix = 64L, n_p = 4L, seed = 12L)
  responses <- data.frame(
    participant_id = sprintf("p%02d", 1:4),
    condition      = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  gcis <- group_ci(sm, responses, by = "condition")
  L <- as.list(gcis)
  expect_named(L, c("A", "B"))
  expect_equal(length(L[[1L]]), nrow(sm))
})

test_that("group_ci output feeds plot_ci_distance_matrix", {
  sm <- make_sig_with_pids(n_pix = 256L, n_p = 8L, seed = 13L)
  responses <- data.frame(
    participant_id = rep(sprintf("p%02d", 1:8), each = 3L),
    condition      = rep(c("A", "B", "C", "D"), each = 6L),
    stringsAsFactors = FALSE
  )
  gcis <- group_ci(sm, responses, by = "condition")
  out <- plot_ci_distance_matrix(gcis,
                                 img_dims    = c(16L, 16L),
                                 show_values = FALSE,
                                 file        = tempfile(fileext = ".pdf"))
  expect_type(out, "list")
  expect_true(is.matrix(out$distance_matrix))
  expect_equal(dim(out$distance_matrix), c(4L, 4L))
})
