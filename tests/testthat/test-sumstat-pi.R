context("SumStat Pi")


test_that("calculation of nuc. div. is correct", {
  seg_sites <- matrix(c(0, 0, 0, 1, 1,
                        0, 1, 1, 0, 0,
                        0, 1, 1, 1, 1), 5, 3)
  expect_equal(calc_nucleotide_div(list(seg_sites), 1:5), 1.6)
  expect_equal(calc_nucleotide_div(list(seg_sites, matrix(0, 5, 0)), 1:5),
               c(1.6, 0.0))
  expect_equal(calc_nucleotide_div(list(seg_sites), 1:3), 4 / 3)
  expect_equal(calc_nucleotide_div(list(seg_sites), c(1, 3, 5)), 2)

  seg_sites <- list(matrix(c(0, 0, 0, 0,
                             0, 0, 1, 1,
                             0, 1, 1, 0,
                             1, 0, 1, 0), 4, 4, TRUE))
  expect_equal(calc_nucleotide_div(seg_sites, 1:4), 2)
})


test_that("nuc. div. statistik works", {
  seg_sites <- matrix(c(0, 0, 0, 1, 1,
                        0, 1, 1, 0, 0,
                        0, 1, 1, 1, 1), 5, 3)
  model <- coal_model(5, 1)
  pi <- sumstat_nucleotide_div("pi", 1)
  expect_equal(pi$calculate(list(seg_sites), NULL, NULL, model), 1.6)

  expect_error(sumstat_nucleotide_div("pi", 1:2))
  expect_error(sumstat_nucleotide_div("pi", 2)$calculate(list(seg_sites),
                                                         NULL, NULL, model))
})
