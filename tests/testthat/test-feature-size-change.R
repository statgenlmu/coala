context("Feature Size Change")

test_that("generating scrm cmd for size changes works", {
  model <- coal_model(15, 1) + feat_size_change(par_range("alpha", 0, 1), 1)
  expect_equal(get_simulator("scrm")$get_cmd(model), "scrm 15 1 -en 0 1 alpha ")
  model <- coal_model(15, 1) + feat_size_change(5, 1)
  expect_equal(get_simulator("scrm")$get_cmd(model), "scrm 15 1 -en 0 1 5 ")
})


test_that("simulating a size change works", {
  scrm <- get_simulator("scrm")
  model_tmp <- coal_model(5, 1) +
    feat_mutation(1) +
    feat_size_change(2, population = 1) +
    sumstat_sfs()

  sum_stats <- scrm$simulate(model_tmp)
  expect_true(is.numeric(sum_stats$sfs))
})
