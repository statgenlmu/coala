context("SumStat iHH")

seg_sites <- matrix(c(1, 0, 0, 0,
                      1, 1, 0, 1,
                      1, 0, 0, 1,
                      1, 0, 0, 1), 4, 4, byrow=TRUE)
attr(seg_sites, 'positions') = c(0.1, 0.2, 0.5, 0.7)
model <- CoalModel(4, 1, 337)


test_that("generation of SNP maps works", {
  stat_iHH <- sumstat_iHH()
  snp_map <- stat_iHH$segsites_to_snp_map(seg_sites, model)
  map <- read.table(snp_map, row.names = 1)
  expect_equal(ncol(map), 4)
  unlink(snp_map)
})


test_that("generation of haplotye file works", {
  stat_iHH <- sumstat_iHH()

  haplotypes <- stat_iHH$segsites_to_haplo(seg_sites)
  haplo <- read.table(haplotypes, row.names = 1)
  expect_equivalent(haplo, as.data.frame(t(seg_sites)))
  unlink(haplotypes)
})


test_that('generation of rehh data works', {
  stat_iHH <- sumstat_iHH()
  rehh_data <- stat_iHH$segsites_to_rehh_data(seg_sites, model)
})


test_that('calculation of iHH works', {
  stat_iHH <- sumstat_iHH()
  iHH <- stat_iHH$calculate(list(seg_sites), NULL, model)
  expect_that(iHH, is_a('list'))
  expect_equal(length(iHH), 1)
  expect_that(iHH[[1]], is_a('matrix'))
  expect_equal(dim(iHH[[1]]), c(2,4))
})

