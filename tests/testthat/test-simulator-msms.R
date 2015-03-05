context('Simulator msms')

test_that("calling msms works", {
  if (!msms_find_jar(FALSE, TRUE)) skip('msms not installed')
  ms.args <- "5 1 -r 10 100 -t 5 -I 2 3 2 1"
  msms.args <- ""
  set.seed(17)
  out_file <- call_msms(ms.args, msms.args)
  set.seed(17)
  out_file_2 <- call_msms(ms.args, msms.args)
  set.seed(20)
  out_file_3 <- call_msms(ms.args, msms.args)
  expect_equal(file.info(out_file_2)$size, file.info(out_file)$size)
  expect_true(file.info(out_file)$size != file.info(out_file_3)$size)
  unlink(c(out_file, out_file_2, out_file_3))
})


test_that("Creation of parameter enviroment works", {
  par_envir <- create_par_env(model_theta_tau(), c(1,5))
  expect_equal(par_envir[['tau']], 1)
  expect_equal(par_envir[['theta']], 5)

  par_envir <- create_par_env(model_theta_tau(), c(1,5), locus = 17)
  expect_equal(par_envir[['locus']], 17)

  par_envir <- create_par_env(model_theta_tau(), c(1,5),
                                  locus = 23, seed = 115)
  expect_equal(par_envir[['locus']], 23)
  expect_equal(par_envir[['seed']], 115)
})


test_that("generating msms options works", {
  dm <- model_theta_tau() +
    feat_selection(par_const(500), par_const(250),
                   population = 1, at_time = 1.7)

  opts <- paste(eval(parse(text = msms_generate_opts_cmd(dm))),
                collapse = " ")
  expect_equal(grep("-SI 1.7 2", opts), 1)
  expect_equal(grep("-SAA 500", opts), 1)
  expect_equal(grep("-SAa 250", opts), 1)
})


test_that("generating text command works", {
  msms <- get_simulator("msms")
  model2 <- model_theta_tau() +
    feat_selection(par_range('s', 1, 10), par_expr(s),
                   population = 1, at_time = 1.7)

  cmd <- msms$get_cmd(model2)
  expect_equal(grep('^msms', cmd), 1)
  expect_equal(grep('-ms 25 10', cmd), 1)
  expect_equal(grep('-SAA s', cmd), 1)
  expect_equal(grep('-SAa s', cmd), 1)
})


test_that("msms_simulate works", {
  if (!msms_find_jar(FALSE, TRUE)) skip('msms not installed')
  msms <- get_simulator("msms")
  dm <- model_theta_tau() +
    feat_selection(par_const(500), par_const(250),
                   population = 1, at_time = 1.7)

  set.seed(6688)
  sum_stats <- msms$simulate(dm, c(1, 5))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)

  set.seed(6688)
  sum_stats2 <- msms$simulate(dm, c(1, 5))
  expect_equal(sum_stats, sum_stats2)
})


test_that("msms_simulate works with inter-locus variation", {
  if (!msms_find_jar(FALSE, TRUE)) skip('msms not installed')
  msms <- get_simulator("msms")

  dm_tmp <- coal_model(5, 2) +
    feat_mutation(par_range('theta', 1, 5), variance = 17) +
    sumstat_seg_sites()
  expect_true(has_inter_locus_var(dm_tmp))

  set.seed(1100)
  sum_stats <- msms$simulate(dm_tmp, c(3))
  expect_is(sum_stats$seg_sites, 'list')
  expect_equal(length(sum_stats$seg_sites), 2)
})


test_that('simulating unphased data works', {
  if (!msms_find_jar(FALSE, TRUE)) skip('msms not installed')
  msms <- get_simulator("msms")
  model <- model_theta_tau() + feat_unphased(2, 1) + sumstat_seg_sites()
  stats <- msms$simulate(model, c(1,5))
  expect_equal(dim(stats$jsfs), c(11, 16))
  expect_equal(nrow(stats$seg_sites[[1]]), 25)

  model <- model_theta_tau() + feat_unphased(3, 2) + sumstat_seg_sites()
  stats <- msms$simulate(model, c(1,5))
  expect_equal(dim(stats$jsfs), c(21, 31))
  expect_equal(nrow(stats$seg_sites[[1]]), 50)
})
