context("SumStat OmegaPrime")

test_that("calculation is correct", {
  ss <- matrix(c(1, 0, 0, 1,
                 1, 1, 0, 0,
                 1, 0, 1, 0,
                 1, 0, 0, 0), 4, 4, byrow = TRUE)

  seg_sites <- list(cbind(ss, ss, ss))
  attr(seg_sites[[1]], "positions") <- rep(c(0.1, 0.2, 0.5, 0.7), 4)
  attr(seg_sites[[1]], "locus") <- rep(c(-1, 0, 1), each = 4)
  expect_equal(calc_omegaprime(seg_sites, 1:4), c(2 / 12))

  ss <- matrix(c(0, 0, 0, 1,
                 0, 0, 1, 0,
                 0, 0, 1, 0,
                 0, 0, 1, 0), 4, 4, byrow = TRUE)
  seg_sites[[2]] <- cbind(ss, ss, ss)
  attr(seg_sites[[2]], "positions") <- rep(c(0.1, 0.2, 0.5, 0.7), 4)
  attr(seg_sites[[2]], "locus") <- rep(c(-1, 0, 1), each = 4)
  expect_equal(calc_omegaprime(seg_sites, 1:4), c(2 / 12, 4 / 12))

  expect_equal(calc_omegaprime(seg_sites, c(1, 3, 4)), c(2 / 12, 4 / 12))
  expect_equal(calc_omegaprime(seg_sites, c(1)), c(0 / 12, 0 / 12))
  expect_error(calc_omegaprime(seg_sites, 1:5))

  seg_sites <- list(matrix(0, 4, 0))
  attr(seg_sites[[1]], "locus") <- numeric()
  attr(seg_sites[[1]], "position") <- numeric()
  expect_true(is.na(calc_omegaprime(seg_sites, 1:4)))
})


test_that("initialzation of statistic works", {
  ss <- matrix(c(1, 0, 0, 1,
                 1, 1, 0, 0,
                 1, 0, 1, 0,
                 1, 0, 0, 0), 4, 4, byrow = TRUE)

  seg_sites <- list(cbind(ss, ss, ss))
  attr(seg_sites[[1]], "positions") <- rep(c(0.1, 0.2, 0.5, 0.7), 4)
  attr(seg_sites[[1]], "locus") <- rep(c(-1, 0, 1), each = 4)

  stat <- sumstat_omegaprime(name = "omega_prime", 1)
  op <- stat$calculate(seg_sites, NULL, NULL, coal_model(4))
  expect_that(op, is_a("numeric"))
  expect_equal(length(op), 1)
  expect_true(op >= 0 & op <= 1)
})


test_that("simulation of omega prime works", {
  set.seed(125)
  model <- model_trios() + sumstat_omegaprime(name = "omega_prime", 1)
  stats <- simulate(model)
  expect_that(stats$omega_prime, is_a("numeric"))
  expect_equal(length(stats$omega_prime), 1)
  expect_true(all(stats$omega_prime >= 0 & stats$omega_prime <= 1))
})
