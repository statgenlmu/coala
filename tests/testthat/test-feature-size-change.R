context("Feature size change")

test_that('Creation of size change features works', {
  expect_equal(feat_size_change(2, 1)$get_rate(), "2")
  expect_equal(feat_size_change(2, 1)$get_population(), 1)

  feat <- feat_size_change(3, 2, 5)
  expect_equal(feat$get_rate(), "3")
  expect_equal(feat$get_population(), 2)
  expect_equal(feat$get_time(), "5")

  expect_error(feat_size_change(3, "A", 5))
  expect_error(feat_size_change(3, 1:2, 5))
})


test_that("generating ms cmd for size changes works", {
  model <- coal_model(15, 1) + feat_size_change(par_range("alpha", 0, 1), 1)
  expect_equal(get_simulator("ms")$get_cmd(model), "ms 15 1 -en 0 1 alpha ")
  model <- coal_model(15, 1) + feat_size_change(5, 1)
  expect_equal(get_simulator("ms")$get_cmd(model), "ms 15 1 -en 0 1 5 ")
})


test_that("simulating a size change works", {
  ms <- get_simulator("ms")
  model_tmp <- coal_model(5, 1) +
    feat_mutation(1) +
    feat_size_change(2, population = 1) +
    sumstat_sfs()

  sum_stats <- ms$simulate(model_tmp)
  expect_true(is.numeric(sum_stats$sfs))
})
