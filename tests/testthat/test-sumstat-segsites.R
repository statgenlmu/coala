context('SumStat SegSites')

test_that('SegSites statistic works', {
  stat <- sumstat_seg_sites('segsites_test')

  seg_sites <- list(matrix(c(1, 0, 0, 0,
                             1, 1, 0, 1,
                             1, 0, 0, 1,
                             1, 0, 0, 1), 4, 4, byrow=TRUE))
  attr(seg_sites[[1]], 'positions') <- c(0.1, 0.2, 0.5, 0.7)

  expect_equal(stat$get_name(), 'segsites_test')
  expect_equal(stat$calculate(seg_sites, NULL, NULL), seg_sites)
})


test_that("SegSites are convert for trios", {
  seg_sites <- list(matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 0,
                             1, 1, 0, 1, 1, 1, 0, 1, 1,
                             1, 0, 0, 1, 1, 0, 0, 1, 1,
                             1, 0, 0, 1, 1, 0, 0, 1, 0), 4, 9, byrow=TRUE))
  attr(seg_sites[[1]], 'positions') <- 1:9 / 10
  llm <- matrix(c(25, 10, 30, 10, 25, 1), 1, 6)

  seg_sites_trio <- seg_sites
  seg_sites_trio[[1]] <- seg_sites_trio[[1]][ , c(1:2, 4:6, 8:9)]
  attr(seg_sites_trio[[1]], "positions") <- c(.4, 0.8,
                                              c(5, 15, 25) / 30,
                                              .2, .6)
  attr(seg_sites_trio[[1]], "locus") <- c(-1, -1, 0, 0, 0, 1, 1)
  expect_equal(conv_for_trios(seg_sites, llm), seg_sites_trio)

  seg_sites <- list(seg_sites[[1]], seg_sites[[1]])
  llm <- rbind(llm, llm)
  seg_sites_trio <- list(seg_sites_trio[[1]], seg_sites_trio[[1]])
  expect_equal(conv_for_trios(seg_sites, llm), seg_sites_trio)
})
