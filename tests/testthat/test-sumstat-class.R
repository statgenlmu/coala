context('SumStat Class')

test_that('sumstat initialization works', {
  ss1 <- sumstat$new('1')
  expect_true(is.sum_stat(ss1))
  expect_error(ss1$calculate())

  ss1 <- sumstat$new('1', 25)
  expect_true(is.sum_stat(ss1))

  expect_false(is.sum_stat(1:3))
})


test_that('Group assignment works', {
  ss1 <- sumstat$new('1')
  expect_equal(ss1$get_group(), 0)
  ss1 <- sumstat$new('1', group = 2)
  expect_equal(ss1$get_group(), 2)
})


test_that('getting the name works', {
  ss1 <- sumstat$new('test')
  expect_equal(ss1$get_name(), 'test')

  ss1 <- sumstat$new('test', group = 2)
  expect_equal(ss1$get_name(), 'test')
})


test_that('Adding sumstats to a model works', {
  model <- coal_model(5:6, 10) + sumstat$new('1')
  expect_equal(get_summary_statistics(model), '1')
  expect_error(model + sumstat$new('1'))

  model <- model + sumstat$new('2')
  expect_equal(get_summary_statistics(model), c('1', '2'))

  model <- model + sumstat$new('3', group=2)
  expect_equal(get_summary_statistics(model, group = 1), c('1', '2'))
  expect_equal(get_summary_statistics(model, group = 2), c('1', '2', '3'))
})


test_that('Getting sumstats groups works', {
  model <- coal_model(5:6, 10) + sumstat$new('1')
  expect_equal(get_sumstat_groups(model), 0)

  model <- model + sumstat$new('2', 15)
  expect_equal(get_sumstat_groups(model), c(0, 15))

  model <- model + sumstat$new('3', 15)
  expect_equal(get_sumstat_groups(model), c(0, 15))
})


test_that('Calculation of sumstats works', {
  Stat_Sum <- R6::R6Class('Stat_Sum', inherit = sumstat,
    public = list(calculate = function(seg_sites, files, model) {
      sapply(seg_sites, sum)
    })
  )
  model <- coal_model(5:6, 10) + Stat_Sum$new('sum')
  stats <- calc_sumstats(list(1:3, 1:5, 7), '', model)
  expect_equal(stats, list(sum=c(6, 15, 7)))

  model <- model + Stat_Sum$new('sum2')
  stats <- calc_sumstats(list(1:3, 1:5, 7), '', model, 1:2)
  expect_equal(stats, list(pars=1:2, sum=c(6, 15, 7), sum2=c(6, 15, 7)))
})
