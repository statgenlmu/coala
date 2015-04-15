context("SumStat iHH")

seg_sites <- matrix(c(1, 0, 0, 0, 1,
                      1, 1, 0, 1, 0,
                      1, 0, 0, 1, 1,
                      1, 0, 0, 1, 0), 4, 5, byrow=TRUE)
attr(seg_sites, 'positions') <- c(0.1, 0.2, 0.5, 0.7, 0.9)
model <- coal_model(4, 1, 337)
pos <- get_snp_positions(list(seg_sites), model, relative = FALSE)[[1]]


test_that("generation of SNP maps works", {
  skip_on_cran()
  stat_ihh <- sumstat_ihh(population = 1)
  snp_map <- stat_ihh$segsites_to_snp_map(seg_sites, pos)
  map <- read.table(snp_map, row.names = 1)
  expect_equal(nrow(map), 5)
  unlink(snp_map)
})


test_that("generation of haplotye file works", {
  skip_on_cran()
  stat_ihh <- sumstat_ihh(population = 1)
  haplotypes <- stat_ihh$segsites_to_haplo(seg_sites, 1:4)
  haplo <- read.table(haplotypes, row.names = 1)
  expect_equivalent(haplo, as.data.frame(seg_sites))
  unlink(haplotypes)
})


test_that('calculation of ihh works', {
  skip_on_cran()
  stat_ihh <- sumstat_ihh()
  ihh <- stat_ihh$calculate(list(seg_sites), NULL, model)
  expect_that(ihh, is_a('list'))
  expect_equal(length(ihh), 1)
  expect_that(ihh[[1]], is_a('matrix'))
  expect_equal(dim(ihh[[1]]), c(5,3))

  stat_ihh <- sumstat_ihh(position = 0.5)
  ihh2 <- stat_ihh$calculate(list(seg_sites), NULL, model)
  expect_that(ihh2, is_a('list'))
  expect_equal(length(ihh2), 1)
  expect_that(ihh2[[1]], is_a('matrix'))
  expect_equal(dim(ihh2[[1]]), c(1,3))
  expect_equal(ihh[[1]][3, , drop=FALSE], ihh2[[1]])

  model <- coal_model(4, 3, 337)
  ihh2 <- stat_ihh$calculate(list(seg_sites, seg_sites, seg_sites), NULL, model)
  expect_that(ihh2, is_a('list'))
  expect_equal(length(ihh2), 3)
  expect_equal(ihh2[[1]], ihh2[[2]])
  expect_equal(ihh2[[1]], ihh2[[3]])
})


test_that('ihh works with trios', {
  skip_on_cran()
  model <- model_trios()
  stats <- simulate(model)
  ihh <- sumstat_ihh(population = 1)
  stat <- ihh$calculate(stats$seg_sites, NULL, model)
  expect_that(stat, is_a('list'))
  expect_equal(length(stat), 1)
})


test_that("ihh works with empty segsites", {
  skip_on_cran()
  model <- model_trios()
  seg_sites <- list(matrix(0, 5, 0))
  ihh <- sumstat_ihh(population = 1)
  attr(seg_sites[[1]], "positions") <- numeric(0)
  stat <- ihh$calculate(seg_sites, NULL, model)
  expect_equivalent(stat, list(matrix(0, 2, 0)))

  seg_sites <- list(matrix(0, 0, 0))
  attr(seg_sites[[1]], "positions") <- numeric(0)
  stat <- ihh$calculate(seg_sites, NULL, model)
  expect_equivalent(stat, list(matrix(0, 2, 0)))
})
