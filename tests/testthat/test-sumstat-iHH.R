context("SumStat iHS")

seg_sites <- matrix(c(1, 0, 0, 0, 1,
                      1, 1, 0, 1, 0,
                      1, 0, 0, 1, 1,
                      1, 0, 0, 1, 0), 4, 5, byrow = TRUE)
attr(seg_sites, "positions") <- c(0.1, 0.2, 0.5, 0.7, 0.9)
model <- coal_model(4, 1, 337)
pos <- get_snp_positions(list(seg_sites), model, relative = FALSE)[[1]]


test_that("generation of rehh data works", {
  skip_if_not_installed("rehh")
  stat_ihh <- sumstat_ihh(population = 1)
  rehh_data <- stat_ihh$create_rehh_data(seg_sites, pos, 1:4)
  expect_equivalent(rehh_data@haplo, seg_sites + 1)
  expect_equal(rehh_data@position, pos)
  expect_equal(rehh_data@snp.name, as.character(1:5))
  expect_equal(rehh_data@nhap, 4)
  expect_equal(rehh_data@nsnp, 5)

  rehh_data <- stat_ihh$create_rehh_data(matrix(0, 4, 0), numeric(), 1:4)
  expect_equal(rehh_data@haplo, matrix(0, 4, 0))

  rehh_data <- stat_ihh$create_rehh_data(seg_sites, pos, numeric())
  expect_equal(rehh_data@haplo, matrix(0, 0, 5))
})


test_that("selection of snps works", {
  stat_ihh <- sumstat_ihh(population = 1, max_snps = 2)
  rehh_data <- stat_ihh$create_rehh_data(seg_sites, pos, 1:4)
  expect_equal(dim(rehh_data@haplo), c(4, 2))
  expect_equal(rehh_data@nsnp, 2)
  expect_equal(rehh_data@nhap, 4)

  stat_ihh <- sumstat_ihh(population = 1, max_snps = 3)
  rehh_data <- stat_ihh$create_rehh_data(seg_sites, pos, 1:4)
  expect_equal(dim(rehh_data@haplo), c(4, 3))
  expect_equal(rehh_data@nsnp, 3)
})


test_that("calculation of ihh works", {
  skip_if_not_installed("rehh")
  stat_ihh <- sumstat_ihh()
  ihh <- stat_ihh$calculate(list(seg_sites), NULL, NULL, model)
  expect_that(ihh, is_a("list"))
  expect_equal(length(ihh), 1)
  expect_that(ihh[[1]], is_a("matrix"))
  expect_equal(dim(ihh[[1]]), c(5, 3))

  stat_ihh <- sumstat_ihh(position = 0.5)
  ihh2 <- stat_ihh$calculate(list(seg_sites), NULL, NULL, model)
  expect_that(ihh2, is_a("list"))
  expect_equal(length(ihh2), 1)
  expect_that(ihh2[[1]], is_a("matrix"))
  expect_equal(dim(ihh2[[1]]), c(1, 3))
  expect_equivalent(ihh[[1]][3, , drop = FALSE], ihh2[[1]])
  expect_equal(rownames(ihh), rownames(ihh2))

  model <- coal_model(4, 3, 337)
  ihh2 <- stat_ihh$calculate(list(seg_sites, seg_sites, seg_sites),
                             NULL, NULL, model)
  expect_that(ihh2, is_a("list"))
  expect_equal(length(ihh2), 3)
  expect_equal(ihh2[[1]], ihh2[[2]])
  expect_equal(ihh2[[1]], ihh2[[3]])
})


test_that("ihh works with trios", {
  skip_if_not_installed("rehh")
  model <- model_trios()
  stats <- simulate(model)
  ihh <- sumstat_ihh(population = 1)
  stat <- ihh$calculate(stats$seg_sites, NULL, NULL, model)
  expect_that(stat, is_a("list"))
  expect_equal(length(stat), 1)
})


test_that("ihh works with empty segsites", {
  skip_if_not_installed("rehh")
  model <- model_trios()
  seg_sites <- list(matrix(0, 5, 0))
  ihh <- sumstat_ihh(population = 1)
  attr(seg_sites[[1]], "positions") <- numeric(0)
  stat <- ihh$calculate(seg_sites, NULL, NULL, model)
  expect_equal(dim(stat[[1]]), c(0, 3))

  seg_sites <- list(matrix(0, 0, 0))
  attr(seg_sites[[1]], "positions") <- numeric(0)
  stat <- ihh$calculate(seg_sites, NULL, NULL, model)
  expect_equal(dim(stat[[1]]), c(0, 3))
})
