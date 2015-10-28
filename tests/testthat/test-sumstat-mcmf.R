context("SumStat MCMF")

test_that("calculation is correct", {
  ss <- matrix(c(1, 0, 0, 1,
                 0, 1, 0, 0,
                 1, 0, 1, 0,
                 1, 0, 0, 0), 4, 4, byrow = TRUE)

  # No trios
  seg_sites <- list(create_segsites(ss, c(0.1, 0.2, 0.5, 0.7), rep(0, 4)))
  expect_equal(calc_mcmf(seg_sites, 1:4, FALSE), .5)
  expect_equal(calc_mcmf(seg_sites, c(1, 3, 4), FALSE), .5)
  expect_equal(calc_mcmf(seg_sites, 2:4, FALSE), 2/3)
  expect_equal(calc_mcmf(seg_sites, 3:4, FALSE), 1)

  # With trios
  seg_sites <- list(create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4)))

  expect_equal(calc_mcmf(seg_sites, 1:4), c(4 / 12))
  expect_equal(calc_mcmf(seg_sites, 2:4), c(4 / 9))
  expect_equal(calc_mcmf(seg_sites, 3:4), c(2 / 3))

  ss <- matrix(c(0, 0, 0, 1,
                 0, 0, 1, 0,
                 0, 0, 1, 0,
                 0, 0, 1, 0), 4, 4, byrow = TRUE)
  seg_sites[[2]] <- create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4))
  expect_equal(calc_mcmf(seg_sites, 1:4), c(c(4 / 12), c(4 / 6)))
  expect_equal(calc_mcmf(seg_sites, 2:4), c(c(4 / 9), NA))
  expect_error(calc_mcmf(seg_sites, 1:5))

  seg_sites <- list(create_segsites(matrix(0, 4, 0)))
  attr(seg_sites[[1]], "locus") <- numeric()
  attr(seg_sites[[1]], "position") <- numeric()
  expect_true(is.na(calc_mcmf(seg_sites, 1:4)))
})


test_that("initialzation of statistic works", {
  ss <- matrix(c(1, 0, 0, 1,
                 1, 1, 0, 0,
                 1, 0, 1, 0,
                 1, 0, 0, 0), 4, 4, byrow = TRUE)

  seg_sites <- list(create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4)))

  stat <- sumstat_mcmf(population = 1)
  op <- stat$calculate(seg_sites, NULL, NULL, coal_model(4))
  expect_that(op, is_a("numeric"))
  expect_equal(length(op), 1)
  expect_true(op >= 0 & op <= 1)
})


test_that("simulation of MCMF works", {
  set.seed(125)
  model <- model_trios() + sumstat_mcmf(name = "mcmf", population = 1)
  stats <- simulate(model)
  expect_that(stats$mcmf, is_a("numeric"))
  expect_equal(length(stats$mcmf), 1)
  expect_true(all(stats$mcmf >= 0 & stats$mcmf <= 1))
})
