context("Simulation Program scrm")

test_that("simulation with scrm works", {
  model <- model_theta_tau()
  sum_stats <- get_simprog('scrm')$sim_func(model, c(1, 5))
  expect_true(is.list(sum_stats))
  expect_equal(length(sum_stats), 2)
  expect_false(is.null(sum_stats$pars))
  expect_false(is.null(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
})

test_that("simulating files works", {
  folder <- tempfile('scrm_test')
  model <- model_theta_tau() + sumstat_file(folder)
  sum_stats <- get_simprog('scrm')$sim_func(model, c(1, 5))
  expect_true(is.list(sum_stats))
  expect_equal(length(sum_stats), 3)
  expect_false(is.null(sum_stats$file))
  expect_true(is.character(sum_stats$file))
  expect_true(file.exists(sum_stats$file))
  unlink(sum_stats$file)
  unlink(folder, recursive = TRUE)
})


test_that('simulating unphased data works', {
  model <- model_theta_tau() + feat_unphased(2, 1) + sumstat_seg_sites()
  stats <- scrm_simulate(model, c(1,5))
  expect_equal(dim(stats$jsfs), c(11, 16))
  expect_equal(nrow(stats$seg_sites[[1]]), 25)

  model <- model_theta_tau() + feat_unphased(3, 2) + sumstat_seg_sites()
  stats <- scrm_simulate(model, c(1,5))
  expect_equal(dim(stats$jsfs), c(21, 31))
  expect_equal(nrow(stats$seg_sites[[1]]), 50)
})
