context("Simulation Program scrm")

test_that("simulation with scrm works", {
  model <- model_theta_tau()
  sum_stats <- getSimProgram('scrm')$sim_func(model, c(1, 5))
  expect_true(is.list(sum_stats))
  expect_equal(length(sum_stats), 2)
  expect_false(is.null(sum_stats$pars))
  expect_false(is.null(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
})

test_that("simulating files works", {
  model <- model_theta_tau() + sumstat_file()
  sum_stats <- getSimProgram('scrm')$sim_func(model, c(1, 5))
  expect_true(is.list(sum_stats))
  expect_equal(length(sum_stats), 3)
  expect_false(is.null(sum_stats$file))
  expect_true(is.character(sum_stats$file))
  expect_true(file.exists(sum_stats$file))
  unlink(sum_stats$file)
})
