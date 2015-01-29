context('SumStat Class')

test_that('SumStat initialization works', {
  ss1 <- SumStat$new('1')
  expect_true(is.sum_stat(ss1))
  expect_equal(ss1$get_name(), '1')
  expect_equal(ss1$get_group(), 0)
  expect_equal(ss1$get_population(), 0)

  ss1 <- SumStat$new('1', 25, 125)
  expect_true(is.sum_stat(ss1))
  expect_equal(ss1$get_name(), '1')
  expect_equal(ss1$get_population(), 25)
  expect_equal(ss1$get_group(), 125)
})


test_that('Adding SumStats to a model works', {
  model <- CoalModel(5:6, 10) + SumStat$new('1')
  expect_equal(get_summary_statistics(model), '1')
  expect_error(model + SumStat$new('1'))

  model <- model + SumStat$new('2')
  expect_equal(get_summary_statistics(model), c('1', '2'))

  model <- model + SumStat$new('3', group=2)
  expect_equal(get_summary_statistics(model, group = 1), c('1', '2'))
  expect_equal(get_summary_statistics(model, group = 2), c('1', '2', '3'))

  #model <- CoalModel(5:7, 10) +
  #  SumStat$new('0', population=0) +
  #  SumStat$new('1', population=1) +
  #  SumStat$new('1', population=3)
  #expect_equal(get_summary_statistics(model, pop = 1), c('0', '1'))
  #expect_equal(get_summary_statistics(model, pop = 2), c('0'))
  #expect_equal(get_summary_statistics(model, pop = 3), c('0', '1'))
})
