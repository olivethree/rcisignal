test_that("rcisignal_diag_result builds with required fields", {
  r <- rcisignal_diag_result("pass", "Label", "All good.")
  expect_s3_class(r, "rcisignal_diag_result")
  expect_equal(r$status, "pass")
  expect_equal(r$label, "Label")
  expect_equal(r$detail, "All good.")
  expect_equal(r$data, list())
})

test_that("rcisignal_diag_result rejects invalid status", {
  expect_error(rcisignal_diag_result("bad", "x", "y"), "one of")
})

test_that("is_rcisignal_diag_result works", {
  expect_true(is_rcisignal_diag_result(rcisignal_diag_result("pass", "x", "y")))
  expect_false(is_rcisignal_diag_result(list(status = "pass")))
  expect_false(is_rcisignal_diag_result(NULL))
})

test_that("print method runs without error for every status", {
  for (s in c("pass", "warn", "fail", "skip")) {
    r <- rcisignal_diag_result(s, "Label", "Detail line.")
    expect_output(print(r), "Label")
  }
})
