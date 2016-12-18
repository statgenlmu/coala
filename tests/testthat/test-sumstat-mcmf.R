context("SumStat MCMF")

ss <- matrix(c(1, 0, 0, 1,
               0, 1, 0, 0,
               1, 0, 1, 0,
               1, 0, 0, 0), 4, 4, byrow = TRUE)
seg_sites <- list(create_segsites(ss, c(0.1, 0.2, 0.5, 0.7)))


test_that("mcmf is correctly calculate for normal loci", {
  expect_equivalent(calc_mcmf(seg_sites, 1:4, FALSE), matrix(.5, 1))
  expect_equivalent(calc_mcmf(seg_sites, c(1, 3, 4), FALSE), matrix(.5, 1))
  expect_equivalent(calc_mcmf(seg_sites, 2:4, FALSE), matrix(2 / 3, 1))
  expect_equivalent(calc_mcmf(seg_sites, 3:4, FALSE), matrix(1, 1))
})


test_that("mcmf is correctly calculation for locus trios", {
  seg_sites <- list(create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4)))

  expect_equivalent(calc_mcmf(seg_sites, 1:4), matrix(c(4 / 12), 1))
  expect_equivalent(calc_mcmf(seg_sites, 2:4), matrix(c(4 / 9), 1))
  expect_equivalent(calc_mcmf(seg_sites, 3:4), matrix(c(2 / 3), 1))

  ss <- matrix(c(0, 0, 0, 1,
                 0, 0, 1, 0,
                 0, 0, 1, 0,
                 0, 0, 1, 0), 4, 4, byrow = TRUE)
  seg_sites[[2]] <- create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4))
  expect_equivalent(calc_mcmf(seg_sites, 1:4), matrix(c(4 / 12, 4 / 6), 2))
  expect_equivalent(calc_mcmf(seg_sites, 2:4), matrix(c(4 / 9, NA), 2))
  expect_error(calc_mcmf(seg_sites, 1:5))

  seg_sites <- list(create_segsites(matrix(0, 4, 0), numeric()))
  expect_true(is.na(calc_mcmf(seg_sites, 1:4)))
})


test_that("mcmf is correctly calculation for diplod loci", {
  expect_equivalent(calc_mcmf(seg_sites, 1:2, FALSE, ploidy = 2), matrix(.75))
  expect_equivalent(calc_mcmf(seg_sites, 1, FALSE, ploidy = 2), matrix(1))
  expect_equivalent(calc_mcmf(seg_sites, 2, FALSE, ploidy = 2), matrix(1))

  seg_sites[[1]] <- create_segsites(rbind(ss, ss), c(0.1, 0.2, 0.5, 0.7))
  expect_equivalent(calc_mcmf(seg_sites, 1:3, FALSE, ploidy = 2), matrix(.75))
  expect_equivalent(calc_mcmf(seg_sites, 1:4, FALSE, ploidy = 2), matrix(.75))
  expect_error(calc_mcmf(seg_sites, 1:5, FALSE, ploidy = 2))
  expect_error(calc_mcmf(seg_sites, 1:4, FALSE, ploidy = 3))
})


test_that("mcmf is correctly calculation for trioplod loci", {
  seg_sites[[1]] <- create_segsites(rbind(ss, ss), c(0.1, 0.2, 0.5, 0.7))
  expect_equivalent(calc_mcmf(seg_sites, 1:2, FALSE, ploidy = 3), matrix(.75))
})


test_that("initialzation of statistic works", {
  stat <- sumstat_mcmf(population = 1)
  expect_equivalent(stat$calculate(seg_sites, NULL, NULL, coal_model(4)), .5)

  seg_sites <- list(create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4)))
  expect_equivalent(stat$calculate(seg_sites, NULL, NULL,
                                   coal_model(4) + locus_trio()), 1 / 3)
})


test_that("mcmf statistics is correct for diploid models", {
  stat <- sumstat_mcmf(population = 1)
  model <- coal_model(2, ploidy = 2)
  expect_equivalent(stat$calculate(seg_sites, NULL, NULL, model), .75)
  expect_equivalent(stat$calculate(seg_sites, NULL, NULL,
                                   model + feat_unphased(2)), .75)
  expect_equivalent(stat$calculate(seg_sites, NULL, NULL,
                                   model + feat_unphased(1)), 1)
})


test_that("mcmf can be calculated for all populations", {
  stat <- sumstat_mcmf(population = "all")
  model1 <- coal_model(c(2, 2))
  model2 <- coal_model(4)
  expect_equal(stat$calculate(seg_sites, NULL, NULL, model1),
               stat$calculate(seg_sites, NULL, NULL, model2))
})


test_that("mcmf can be calculated for multiple populations", {
  model <- coal_model(c(5, 5, 1), 1, ploidy = 2) +
    sumstat_mcmf("mcmf1", population = 1) +
    sumstat_mcmf("mcmf2", population = 2) +
    feat_mutation(10) +
    feat_migration(1, symmetric = TRUE)
  stats <- simulate(model)
  expect_is(stats$mcmf1, "numeric")
  expect_is(stats$mcmf2, "numeric")

  model <- model + feat_unphased(1)
  stats <- simulate(model)
  expect_is(stats$mcmf1, "numeric")
  expect_is(stats$mcmf2, "numeric")
})


test_that("mcmf can be calulated in the expanded version", {
  model <- coal_model(4) + locus_averaged(1, 4)

  seg_sites <- list(create_segsites(cbind(ss, ss, ss),
                                    rep(c(0.1, 0.2, 0.5, 0.7), 3),
                                    rep(c(-1, 0, 1), each = 4)))

  stat <- sumstat_mcmf(population = 1, expand_mcmf = TRUE, type_expand = 1)
  expect_equivalent(stat$calculate(seg_sites, NULL, NULL, model), matrix(0.5))

  stat <- sumstat_mcmf(population = 1, expand_mcmf = TRUE, type_expand = 2)
  expect_equivalent(stat$calculate(seg_sites, NULL, NULL, model),
                    matrix(c(0.5, 0.25), 1))

  stat <- sumstat_mcmf(population = 1, expand_mcmf = TRUE, type_expand = 3)
  expect_equivalent(stat$calculate(seg_sites, NULL, NULL, model),
                    matrix(c(0.5, 0.25, 1), 1))

  stat <- sumstat_mcmf(population = 1, expand_mcmf = TRUE, type_expand = 3)
  expect_equivalent(stat$calculate(list(seg_sites[[1]], seg_sites[[1]]),
                                   NULL, NULL, model),
                    matrix(c(0.5, 0.25, 1), 2, 3, byrow = TRUE))

  model <- coal_model(4) + locus_trio()
  stat <- sumstat_mcmf(population = 1, expand_mcmf = TRUE, type_expand = 3)
  expect_equivalent(stat$calculate(seg_sites, NULL, NULL, model),
                    matrix(c(1 / 3, 0.25, 0.004), 1))
})
