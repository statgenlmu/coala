context("SumStat FourGamete")


test_that("calc_four_gamete_stat works", {
  seg_sites <- list(matrix(c(1, 1, 0, 0, 0,
                             1, 0, 1, 0, 1,
                             1, 1, 0, 1, 0,
                             0, 1, 1, 0, 1,
                             0, 0, 0, 0, 1), 5))

  attr(seg_sites[[1]], "positions") <- c(0.1, 0.12, 0.5, 0.51, 0.61)
  locus_length <- matrix(c(0, 0, 100, 0, 0, 1), 1, 6)
  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[1, ], c(mid_near = .5, mid_far = .5, outer = NaN,
                                      between = NaN, mid = .5,
                                      perc_polym = 0.05))

  seg_sites[[2]] <- seg_sites[[1]]
  locus_length <- rbind(locus_length, matrix(c(0, 0, 50, 0, 0, 1), 1, 6))
  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[1, ], c(mid_near = .5, mid_far = .5, outer = NaN,
                                      between = NaN, mid = .5,
                                      perc_polym = 0.05))
  expect_equal(fpc_violations[2, ], c(mid_near = .5, mid_far = .5, outer = NaN,
                                      between = NaN, mid = .5,
                                      perc_polym = 0.10))

  attr(seg_sites[[2]], "positions")[4:5] <- c(0.7, 0.75)
  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[2, ], c(mid_near = 1, mid_far = 0.4, outer = NaN,
                                      between = NaN, mid = 0.5,
                                      perc_polym = 0.10))

  attr(seg_sites[[1]], "positions") <- 1:5 / 5
  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[1, ], c(mid_near = NaN, mid_far = 0.5,
                                      outer = NaN, between = NaN, mid = 0.5,
                                      perc_polym = 0.05))
  attr(seg_sites[[1]], "positions") <- 1:5 / 50
  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[1, ], c(mid_near = 0.5, mid_far = NaN,
                                      outer = NaN, between = NaN, mid = 0.5,
                                      perc_polym = 0.05))

  seg_sites[[2]][1, ] <- 1
  seg_sites[[2]][-1, ] <- 0
  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[2, ], c(mid_near = NA, mid_far = NA, outer = NA,
                                      between = NA, mid = NA,
                                      perc_polym = 0.10))


  seg_sites[[3]] <- matrix(0, 5, 0)
  attr(seg_sites[[3]], "positions") <- numeric(0)
  locus_length <- rbind(locus_length, c(0, 0, 50, 0, 0, 1))
  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[3, ], c(mid_near = NaN, mid_far = NaN,
                                      outer = NaN, between = NaN, mid = NaN,
                                      perc_polym = 0))


  # With locus-trios
  seg_sites[[4]] <- matrix(c(1, 1, 0, 0, 0,
                             1, 0, 1, 0, 1,
                             1, 1, 1, 1, 0,
                             0, 1, 1, 0, 1,
                             0, 0, 0, 0, 1), 5)
  attr(seg_sites[[4]], "positions") <- c(0.15, 0.55, 0.05, 0.08, 0.30)
  attr(seg_sites[[4]], "locus") <- c(-1, -1, 0, 0, 0)
  locus_length <- rbind(locus_length, c(10, 5, 6, 5, 10, 1))

  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[4, ], c(mid_near = NaN, mid_far = NaN, outer = 1,
                                      between = 1, mid = NaN, perc_polym = 0.5))

  attr(seg_sites[[4]], "locus") <- c(0, 0, 1, 1, 1)
  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[4, ], c(mid_near = NaN, mid_far = 1, outer = NaN,
                                      between = 1, mid = 1, perc_polym = 1 / 3))

  attr(seg_sites[[4]], "locus") <- c(-1, -1, 0, 0, 1)
  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[4, ], c(mid_near = NaN, mid_far = NaN, outer = 1,
                                      between = 1, mid = NaN,
                                      perc_polym = 1 / 3))

  attr(seg_sites[[4]], "locus") <- c(-1, -1, 0, 1, 1)
  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[4, ], c(mid_near = NaN, mid_far = NaN, outer = 1,
                                      between = NaN, mid = NaN,
                                      perc_polym = 1 / 6))

  attr(seg_sites[[4]], "locus") <- c(-1, -1, -1, 1, 1)
  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[4, ], c(mid_near = NaN, mid_far = NaN, outer = 1,
                                      between = NaN, mid = NaN, perc_polym = 0))

  # Not all individuals
  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:2, locus_length)
  expect_equal(fpc_violations[4, ], c(mid_near = NaN, mid_far = NaN, outer = NaN,
                                      between = NaN, mid = NaN,
                                      perc_polym = 0))

  fpc_violations <- calc_four_gamete_stat(seg_sites, 1:4, locus_length)
  expect_equal(fpc_violations[4, ], c(mid_near = NaN, mid_far = NaN, outer = 1,
                                      between = NaN, mid = NaN,
                                      perc_polym = 0))
})


test_that("Simulation the statistic works", {
  model <- model_theta_tau() + sumstat_four_gamete()
  stats <- simulate(model, pars = c(1, 5))
  expect_that(stats$four_gamete, is_a("matrix"))
  expect_equal(nrow(stats$four_gamete), get_locus_number(model))
})
