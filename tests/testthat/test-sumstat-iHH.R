context("SumStat iHH")

seg_sites <- matrix(c(1, 0, 0, 0, 1,
                      1, 1, 0, 1, 0,
                      1, 0, 0, 1, 1,
                      1, 0, 0, 1, 0), 4, 5, byrow=TRUE)
attr(seg_sites, 'positions') = c(0.1, 0.2, 0.5, 0.7, 0.9)
model <- CoalModel(4, 1, 337)
pos <- get_snp_positions(list(seg_sites), model, relative = FALSE)[[1]]


test_that("generation of SNP maps works", {
  stat_iHH <- sumstat_iHH(population = 1)
  snp_map <- stat_iHH$segsites_to_snp_map(seg_sites, pos)
  map <- read.table(snp_map, row.names = 1)
  expect_equal(nrow(map), 5)
  unlink(snp_map)
})


test_that("generation of haplotye file works", {
  stat_iHH <- sumstat_iHH(population = 1)
  haplotypes <- stat_iHH$segsites_to_haplo(seg_sites, 1:4)
  haplo <- read.table(haplotypes, row.names = 1)
  expect_equivalent(haplo, as.data.frame(seg_sites))
  unlink(haplotypes)
})


test_that('calculation of iHH works', {
  stat_iHH <- sumstat_iHH()
  iHH <- stat_iHH$calculate(list(seg_sites), NULL, model)
  expect_that(iHH, is_a('list'))
  expect_equal(length(iHH), 1)
  expect_that(iHH[[1]], is_a('matrix'))
  expect_equal(dim(iHH[[1]]), c(2,5))

  stat_iHH <- sumstat_iHH(position = 0.5)
  iHH2 <- stat_iHH$calculate(list(seg_sites), NULL, model)
  expect_that(iHH2, is_a('list'))
  expect_equal(length(iHH2), 1)
  expect_that(iHH2[[1]], is_a('matrix'))
  expect_equal(dim(iHH2[[1]]), c(2,1))
  expect_equal(iHH[[1]][,3,drop=FALSE], iHH2[[1]])
})


test_that('iHH works with trios', {
  model <- model_trios()
  stats <- simulate(model, pars=c(1,5))
  ihh <- sumstat_iHH(population = 2)
  ihh$calculate(stats$seg_sites, NULL, model)
})
