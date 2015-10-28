context("SumStat Omega")

test_that("Initialization of Omega works", {
  if (!has_omega()) skip("OmegaPlus not found")
  op <- sumstat_omega(name = "op", min_win = 12, max_win = 112, grid = 15)
  expect_true(is.sum_stat(op))
  expect_equal(op$get_name(), "op")
  expect_equal(op$get_min_win(), 12)
  expect_equal(op$get_max_win(), 112)
  expect_equal(op$get_grid(), 15)

  expect_error(sumstat_omega(binary = tempfile("op")))
})


test_that("report files are parsed correctly", {
  if (!has_omega()) skip("OmegaPlus not found")
  tmp_dir <- tempfile("op_parse_test")
  dir.create(tmp_dir)

  cat("//1
1.11	1.000000
2.22	0.200000
3.33	0.030000

//2
1.23	4.000000
2.34	0.500000
3.45	0.060000
", file = file.path(tmp_dir, "OmegaPlus_Report.op"))

  op <- sumstat_omega(name = "op", min_win = 12, max_win = 112, grid = 5)
  expect_equal(op$parse_report(tmp_dir, n_grid = 3, start_locus = 1),
               data.frame(locus = c(1, 1, 1, 2, 2, 2),
                          pos = c(1.11, 2.22, 3.33, 1.23, 2.34, 3.45),
                          omega = c(1, .2, .03, 4, .5, .06)))
  unlink(tmp_dir, recursive = TRUE)
})


test_that("Omega can be calculate", {
  if (!has_omega()) skip("OmegaPlus not found")
  model <- coal_model(10, 2) +
    feat_mutation(5) +
    sumstat_omega("op", grid = 10)
  stat <- simulate(model)
  expect_false(is.null(stat$op))
  expect_equal(dim(stat$op), c(20, 3))
})


test_that("Omega checks that the number of grid points is valid", {
  if (!has_omega()) skip("OmegaPlus not found")
  expect_error(coal_model(10, 2, 100) +
                 feat_mutation(5) +
                 sumstat_omega("op", grid = 1000))
})


test_that("Omega works if there are few SNPs", {
  if (!has_omega()) skip("OmegaPlus not found")
  model <- coal_model(10, 2, 1000) +
    feat_mutation(2, fixed_number = TRUE) +
    sumstat_omega("op", grid = 10)
  simulate(model)
})


test_that("OmegaPrime rejects trio loci", {
  if (!has_omega()) skip("OmegaPlus not found")
  expect_error(coal_model(10) +
    feat_mutation(5) +
    sumstat_omega("op") +
    locus_trio(number = 2))
})
