context("SumStat XP-CLR")

test_seg_sites <- create_segsites(matrix(c(1, 0, 0, 0, 1,
                                           1, 1, 0, 1, 0,
                                           1, 0, 1, 1, 1,
                                           0, 0, 0, 1, 0), 4, 5, byrow = TRUE),
                                  c(0.1, 0.2, 0.5, 0.7, 0.9))

test_model <- coal_model(c(2, 2), 2, 337)


test_that("it creates geno files", {
  if (!has_xp_clr()) skip("XPCLR not found")
  xp_clr <- sumstat_xp_clr("xp-clr", 1, 2)
  geno_file <- xp_clr$create_geno_file(list(test_seg_sites, test_seg_sites), test_model, 1)
  expect_true(file.exists(geno_file))

  geno_data <- as.matrix(read.table(geno_file))
  expect_equivalent(geno_data, matrix(c(1, 0, 0, 0, 1, 1, 0, 0, 0, 1,
                                        1, 1, 0, 1, 0, 1, 1, 0, 1, 0), 10, 2))
  unlink(geno_file)

  geno_file <- xp_clr$create_geno_file(list(test_seg_sites), test_model, 2)
  expect_true(file.exists(geno_file))

  geno_data <- as.matrix(read.table(geno_file))
  expect_equivalent(geno_data, matrix(c(1, 0, 1, 1, 1,
                                        0, 0, 0, 1, 0), 5, 2))
  unlink(geno_file)
})

