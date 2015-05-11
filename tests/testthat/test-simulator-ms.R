context('Simulator ms')

test_that("the ms sim program exists", {
  expect_that(get_simulator("ms"), is_a("Simulator"))
})


test_that("msSimFunc is working", {
  ms <- get_simulator("ms")
  model_tt <- model_theta_tau()
  set.seed(789)
  sum_stats <- ms$simulate(model_tt, c(tau = 1, theta = 10))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
})


test_that("msSimFunc works with inter-locus variation", {
  ms <- get_simulator("ms")
  model_tmp <- coal_model(5:6, 2) +
    feat_mutation(par_variation(2, 2)) +
    feat_pop_merge(par_const(.5), 2, 1) +
    sumstat_jsfs()

  set.seed(117)
  sum_stats <- ms$simulate(model_tmp)
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)

  set.seed(117)
  sum_stats2 <- ms$simulate(model_tmp)
  expect_equal(sum_stats$jsfs, sum_stats2$jsfs)
})


test_that('simulating unphased data works', {
  ms <- get_simulator("ms")
  model <- model_theta_tau() + feat_unphased(2, 1) + sumstat_seg_sites()
  stats <- ms$simulate(model, c(tau = 1, theta = 5))
  expect_equal(dim(stats$jsfs), c(11, 16))
  expect_equal(nrow(stats$seg_sites[[1]]), 25)

  model <- model_theta_tau() + feat_unphased(3, 2) + sumstat_seg_sites()
  stats <- ms$simulate(model, c(tau = 1, theta = 5))
  expect_equal(dim(stats$jsfs), c(21, 31))
  expect_equal(nrow(stats$seg_sites[[1]]), 50)
})


test_that("ms can simulate locus trios", {
  stat <- get_simulator("ms")$simulate(model_trios())
  expect_that(attr(stat$seg_sites[[1]], "locus"), is_a("numeric"))
  expect_true(all(attr(stat$seg_sites[[1]], "locus") %in% -1:1))
  expect_true(all(attr(stat$seg_sites[[1]], "positions") >= 0))
  expect_true(all(attr(stat$seg_sites[[1]], "positions") <= 1))
})


test_that("ms works with scientific notation", {
  model <- coal_model(5, 1, 1e8) + feat_recombination(1)
  opts <- ms_generate_opts(model, numeric(0), 1)
  expect_true("100000000" %in% opts)
  expect_false("1e8" %in% opts)

  model <- coal_model(5, 1, 1000) + feat_recombination(1e8)
  opts <- ms_generate_opts(model, numeric(0), 1)
  expect_true("100000000" %in% opts)
  expect_false("1e+08" %in% opts)

  model <- coal_model(5, 1, 1000) + feat_mutation(1e8)
  opts <- ms_generate_opts(model, numeric(0), 1)
  expect_true("100000000" %in% opts)
  expect_false("1e+08" %in% opts)
})
