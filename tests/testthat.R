if (require("testthat")) {
  test_check("coala")
} else {
  warning("testthat not available. Skipping unittests!")
}
