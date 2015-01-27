context("Simulation Program scrm")

test_that("simulation with scrm works", {
  dm_tt <- model_theta_tau()
  sum_stats <- getSimProgram('scrm')$sim_func(dm_tt, c(1, 5))
  expect_true(is.list(sum_stats))
  expect_equal(length(sum_stats), 2)
  expect_false(is.null(sum_stats$pars))
  expect_false(is.null(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
})

test_that("simulating files works", {
  dm_tt <- model_theta_tau()
  dm <- dm.addSummaryStatistic(dm_tt, 'file')
  sum_stats <- getSimProgram('scrm')$sim_func(dm, c(1, 5))
  expect_true(is.list(sum_stats))
  expect_equal(length(sum_stats), 3)
  expect_false(is.null(sum_stats$file))
  expect_true(is.character(sum_stats$file))
  expect_true(file.exists(sum_stats$file))
  unlink(sum_stats$file)
})
