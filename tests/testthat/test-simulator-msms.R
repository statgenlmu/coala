context("Simulator msms")

test_that("msms can be added manually", {
  if (!has_msms()) skip("msms not installed")
  msms_jar <- get_simulator("msms")$get_info()["jar"]
  java <- get_simulator("msms")$get_info()["java"]
  activate_msms(msms_jar, java, 199)
  expect_equal(get_simulator("msms")$get_priority(), 199)
  expect_error(use_msms(tempfile("not-existant"), tempfile("not-existant")))
})


test_that("msms supports size changes in one pop models", {
  if (!has_msms()) skip("msms not installed")

  model <- coal_model(40, 1) +
    feat_mutation(1) +
    feat_size_change(0.1, population = 1, time = 0.1) +
    feat_selection(1000, time = 0.01) +
    sumstat_sfs()

  stat <- simulate(model)
  expect_that(stat$sfs, is_a("numeric"))
})


test_that("msms supports growth in one pop models", {
  if (!has_msms()) skip("msms not installed")

  model <- coal_model(40, 1) +
    feat_mutation(1) +
    feat_growth(0.1, population = 1, time = 0.1) +
    feat_selection(1000, time = 0.01) +
    sumstat_sfs()

  stat <- simulate(model)
  expect_that(stat$sfs, is_a("numeric"))
})
