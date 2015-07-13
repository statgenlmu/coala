context("SumStat SFS")

test_that("calculation of the SFS is correct", {
  seg.sites <- list(matrix(c(1, 0, 0, 0,
                             1, 1, 0, 1,
                             0, 0, 1, 1,
                             1, 1, 0, 1), 4, 4, byrow = TRUE))
  attr(seg.sites[[1]], "positions") <- c(0.1, 0.2, 0.5, 0.7)

  model <- coal_model(4, 1)
  stat_sfs_all <- sumstat_sfs(population = "all")
  stat_sfs_1 <- sumstat_sfs(population = 1)
  stat_sfs_2 <- sumstat_sfs(population = 2)

  expect_equal(stat_sfs_1$calculate(seg.sites, NULL, NULL, model), c(1, 1, 2))
  expect_equal(stat_sfs_all$calculate(seg.sites, NULL, NULL, model), c(1, 1, 2))

  model <- coal_model(c(2, 2), 1)
  expect_equal(stat_sfs_all$calculate(seg.sites, NULL, NULL, model), c(1, 1, 2))
  expect_equal(stat_sfs_1$calculate(seg.sites, NULL, NULL, model), 2)
  expect_equal(stat_sfs_2$calculate(seg.sites, NULL, NULL, model), 3)

  seg.sites[[2]] <- matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),  4, 3)
  attr(seg.sites[[2]], "positions") <- c(0.1, 0.5, 0.7)
  expect_equal(stat_sfs_all$calculate(seg.sites, NULL, NULL, model), c(1, 1, 2))
  expect_equal(stat_sfs_1$calculate(seg.sites, NULL, NULL, model), c(2))
  expect_equal(stat_sfs_2$calculate(seg.sites, NULL, NULL, model), c(3))

  seg.sites[[3]] <- matrix(numeric(), 4, 0)
  attr(seg.sites[[3]], "positions") <- c()
  expect_equal(stat_sfs_all$calculate(seg.sites, NULL, NULL, model), c(1, 1, 2))
  expect_equal(stat_sfs_1$calculate(seg.sites, NULL, NULL, model), 2)
  expect_equal(stat_sfs_2$calculate(seg.sites, NULL, NULL, model), 3)
})


test_that("calculation of sfs works with trios", {
  ss <- matrix(c(1, 0, 0, 0,
                 1, 1, 0, 1,
                 1, 0, 0, 1,
                 1, 0, 0, 1), 4, 4, byrow = TRUE)

  seg.sites <- list(cbind(ss, ss, ss))
  attr(seg.sites[[1]], "positions") <- rep(c(0.1, 0.2, 0.5, 0.7), 4)
  attr(seg.sites[[1]], "locus") <- rep(c(-1, 0, 1), each = 4)

  model <- coal_model(4, 1)
  stat_sfs_all <- sumstat_sfs(population = "all")
  expect_equal(stat_sfs_all$calculate(seg.sites, NULL, NULL, model), c(1, 0, 1))
})
