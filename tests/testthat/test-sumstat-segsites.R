context('SumStat SegSites')

test_that('SegSites statistic works', {
  stat <- sumstat_seg_sites('segsites_test')

  seg_sites <- list(matrix(c(1, 0, 0, 0,
                             1, 1, 0, 1,
                             1, 0, 0, 1,
                             1, 0, 0, 1), 4, 4, byrow=TRUE))
  attr(seg_sites[[1]], 'positions') = c(0.1, 0.2, 0.5, 0.7)

  expect_equal(stat$get_name(), 'segsites_test')
  expect_equal(stat$calculate(seg_sites, NULL, NULL), seg_sites)
})
