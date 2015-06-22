context('Simulator msms')

test_that("calling msms works", {
  if (!has_msms()) skip('msms not installed')
  msms <- get_simulator("msms")
  msms.args <- "5 1 -r 10 100 -t 5 -I 2 3 2 1"
  set.seed(17)
  out_file <- msms$call_msms(msms.args)
  set.seed(17)
  out_file_2 <- msms$call_msms(msms.args)
  set.seed(20)
  out_file_3 <- msms$call_msms(msms.args)
  expect_equal(file.info(out_file_2)$size, file.info(out_file)$size)
  expect_true(file.info(out_file)$size != file.info(out_file_3)$size)
  unlink(c(out_file, out_file_2, out_file_3))
})


test_that("generating msms options works", {
  if (!has_msms()) skip('msms not installed')
  msms <- get_simulator("msms")
  model <- coal_model(10, 2) + feat_mutation(5)
  expect_equal(msms$get_cmd(model), "msms 10 2 -t 5 ")
})


test_that("msms_simulate works", {
  if (!has_msms()) skip('msms not installed')
  msms <- get_simulator("msms")
  model <- model_theta_tau()

  set.seed(6688)
  sum_stats <- msms$simulate(model, c(tau = 1, theta = 5))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)

  set.seed(6688)
  sum_stats2 <- msms$simulate(model, c(tau = 1, theta = 5))
  expect_equal(sum_stats, sum_stats2)

  # With recombination
  model <- model + feat_recombination(1)
  sum_stats <- msms$simulate(model, c(tau = 1, theta = 5))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
})


test_that("msms_simulate works with inter-locus variation", {
  if (!has_msms()) skip('msms not installed')
  msms <- get_simulator("msms")

  model_tmp <- coal_model(5, 2) +
    feat_mutation(par_variation(par_range('theta', 1, 5), 17)) +
    sumstat_seg_sites()
  expect_true(has_variation(model_tmp))

  sum_stats <- msms$simulate(model_tmp, c(theta = 3))
  expect_is(sum_stats$seg_sites, 'list')
  expect_equal(length(sum_stats$seg_sites), 2)
})


test_that('simulating unphased data works', {
  if (!has_msms()) skip('msms not installed')
  msms <- get_simulator("msms")
  model <- model_theta_tau() + feat_unphased(2, 1) + sumstat_seg_sites()
  stats <- msms$simulate(model, c(tau = 1, theta = 5))
  expect_equal(dim(stats$jsfs), c(11, 16))
  expect_equal(nrow(stats$seg_sites[[1]]), 25)

  model <- model_theta_tau() + feat_unphased(3, 2) + sumstat_seg_sites()
  stats <- msms$simulate(model, c(tau = 1, theta = 5))
  expect_equal(dim(stats$jsfs), c(21, 31))
  expect_equal(nrow(stats$seg_sites[[1]]), 50)
})


test_that("msms can simulate locus trios", {
  if (!has_msms()) skip('msms not installed')
  stat <- get_simulator("msms")$simulate(model_trios())
  expect_that(attr(stat$seg_sites[[1]], "locus"), is_a("numeric"))
  expect_true(all(attr(stat$seg_sites[[1]], "locus") %in% -1:1))
  expect_true(all(attr(stat$seg_sites[[1]], "positions") >= 0))
  expect_true(all(attr(stat$seg_sites[[1]], "positions") <= 1))
})


test_that("ms can added manually", {
  if (!has_msms()) skip("msms not installed")
  msms_jar <- get_simulator("msms")$get_info()["jar"]
  java <- get_simulator("msms")$get_info()["java"]
  activate_msms(msms_jar, java, 199)
  expect_equal(get_simulator("msms")$get_priority(), 199)
  expect_error(use_msms(tempfile("not-existant"), tempfile("not-existant")))
})
