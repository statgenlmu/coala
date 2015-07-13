context("SumStat Class")

test_that("Sumstat initialization works", {
  ss1 <- sumstat_class$new("1")
  expect_true(is.sum_stat(ss1))
  expect_error(ss1$calculate())

  ss1 <- sumstat_class$new("1")
  expect_true(is.sum_stat(ss1))

  expect_false(is.sum_stat(1:3))
})


test_that("getting the name works", {
  ss1 <- sumstat_class$new("test")
  expect_equal(ss1$get_name(), "test")
})


test_that("Adding Sumstats to a model works", {
  model <- coal_model(5:6, 10) + sumstat_class$new("1")
  expect_equal(get_summary_statistics(model)[[1]]$get_name(), "1")
  expect_error(model + sumstat_class$new("1"))

  model <- model + sumstat_class$new("2")
  expect_equal(names(get_summary_statistics(model)), c("1", "2"))
})


test_that("Calculation of Sumstats works", {
  Stat_Sum <- R6::R6Class("Stat_Sum", inherit = sumstat_class,
    public = list(calculate = function(seg_sites, trees, files, model) {
      sapply(seg_sites, sum)
    })
  )
  model <- coal_model(5:6, 10) + Stat_Sum$new("sum")
  stats <- calc_sumstats(list(1:3, 1:5, 7), NULL, "", model, 1,
                         1:3, get_simulator("scrm"))
  expect_equal(stats$sum, c(6, 15, 7))

  model <- model + Stat_Sum$new("sum2")
  stats <- calc_sumstats(list(1:3, 1:5, 7), NULL, "", model, 1:2,
                         1:3, get_simulator("scrm"))
  expect_equal(stats$pars, 1:2)
  expect_equal(stats$sum,  c(6, 15, 7))
  expect_equal(stats$sum2, c(6, 15, 7))
  expect_equal(stats$cmds, 1:3)
})
