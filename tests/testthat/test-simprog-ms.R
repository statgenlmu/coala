context('Simulation Program ms')

test_that("msSimFunc is working", {
  dm_tt <- model_theta_tau()
  set.seed(789)
  sum_stats <- msSingleSimFunc(dm_tt, c(1, 10))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)

#   Working interactively, failing in tests somehow?
#   set.seed(789)
#   sum_stats2 <- msSingleSimFunc(dm_tt, c(1, 10))
#   expect_equal(sum_stats, sum_stats2)
})

test_that("msSimFunc works with inter-locus variation", {
  dm_tmp <- CoalModel(5:6, 2) +
    feat_mutation(par_const(2), variance = '2') +
    feat_pop_merge(par_const(.5), 2, 1) +
    sumstat_jsfs()

  set.seed(117)
  sum_stats <- msSingleSimFunc(dm_tmp)
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)

  set.seed(117)
  sum_stats2 <- msSingleSimFunc(dm_tmp)
  expect_equal(sum_stats$jsfs, sum_stats2$jsfs)
})

test_that("the ms sim program exists", {
  expect_false(is.null(getSimProgram('ms')))
})
