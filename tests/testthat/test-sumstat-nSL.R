context("SumStat nSL")

test_that('calculation of nSL works', {
  skip_on_cran()
  seg_sites <- matrix(c(1, 0, 0, 0, 1,
                        1, 1, 0, 1, 0,
                        1, 0, 0, 1, 1,
                        1, 0, 0, 1, 0), 4, 5, byrow = TRUE)
  attr(seg_sites, 'positions') <- c(0.1, 0.2, 0.5, 0.7, 0.9)
  model <- coal_model(4, 1, 337)
  pos <- get_snp_positions(list(seg_sites), model, relative = FALSE)[[1]]

  stat_nsl <- sumstat_nSL()
  nsl <- stat_nsl$calculate(list(seg_sites), NULL, model)
  expect_that(nsl, is_a('list'))
  expect_equal(length(nsl), 1)
  expect_that(nsl[[1]], is_a('numeric'))
  expect_equal(length(nsl[[1]]), 5)

  stat_nsl <- sumstat_nSL(position = 0.5)
  nsl2 <- stat_nsl$calculate(list(seg_sites), NULL, model)
  expect_that(nsl2, is_a('list'))
  expect_equal(length(nsl2), 1)
  expect_that(nsl2[[1]], is_a('numeric'))
  expect_equal(length(nsl2[[1]]), 1)
  expect_equal(nsl[[1]][3], nsl2[[1]])

  model <- coal_model(4, 3, 337)
  nsl <- stat_nsl$calculate(list(seg_sites, seg_sites, seg_sites), NULL, model)
  expect_that(nsl, is_a('list'))
  expect_equal(length(nsl), 3)
  expect_equal(nsl[[1]], nsl[[2]])
  expect_equal(nsl[[1]], nsl[[3]])
})


test_that("nSL works with empty segsites", {
  skip_on_cran()
  model <- model_trios()
  seg_sites <- list(matrix(0, 5, 0))
  nsl <- sumstat_nSL(population = 1)
  attr(seg_sites[[1]], "positions") <- numeric(0)
  stat <- nsl$calculate(seg_sites, NULL, model)
  expect_equal(length(stat[[1]]), 0)

  seg_sites <- list(matrix(0, 0, 0))
  attr(seg_sites[[1]], "positions") <- numeric(0)
  stat <- nsl$calculate(seg_sites, NULL, model)
  expect_equal(length(stat[[1]]), 0)
})
