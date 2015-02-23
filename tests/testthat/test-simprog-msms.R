context('Simulation Program msms')

test_that("calling msms works", {
  if (!checkForMsms(FALSE, TRUE)) skip('msms not installed')
  ms.args <- "5 1 -r 10 100 -t 5 -I 2 3 2 1"
  msms.args <- ""
  set.seed(17)
  out.file <- callMsms(ms.args, msms.args)
  set.seed(17)
  out.file.2 <- callMsms(ms.args, msms.args)
  set.seed(20)
  out.file.3 <- callMsms(ms.args, msms.args)
  expect_equal(file.info(out.file.2)$size, file.info(out.file)$size)
  expect_true(file.info(out.file)$size != file.info(out.file.3)$size)
  unlink(c(out.file, out.file.2, out.file.3))
})


test_that("Creation of parameter enviroment works", {
  par_envir <- createParameterEnv(model_theta_tau(), c(1,5))
  expect_equal(par_envir[['tau']], 1)
  expect_equal(par_envir[['theta']], 5)

  par_envir <- createParameterEnv(model_theta_tau(), c(1,5), locus = 17)
  expect_equal(par_envir[['locus']], 17)

  par_envir <- createParameterEnv(model_theta_tau(), c(1,5), locus = 23, seed = 115)
  expect_equal(par_envir[['locus']], 23)
  expect_equal(par_envir[['seed']], 115)
})


test_that("generating msms options works", {
  dm <- model_theta_tau() +
    feat_selection(par_const(500), par_const(250), population = 1, at_time = 1.7)

  opts <- paste(eval(parse(text = generateMsmsOptionsCommand(dm))),
                collapse = " ")
  expect_equal(grep("-SI 1.7 2", opts), 1)
  expect_equal(grep("-SAA 500", opts), 1)
  expect_equal(grep("-SAa 250", opts), 1)
})


test_that("generating text command works", {
  model2 <- model_theta_tau() +
    feat_selection(par_range('s', 1, 10), par_expr(s),
                   population = 1, at_time = 1.7)

  cmd <- msms_get_command(model2)
  expect_equal(grep('^msms', cmd), 1)
  expect_equal(grep('-ms 25 10', cmd), 1)
  expect_equal(grep('-SAA s', cmd), 1)
  expect_equal(grep('-SAa s', cmd), 1)
})


test_that("msmsSimFunc works", {
  if (!checkForMsms(FALSE, TRUE)) skip('msms not installed')
  dm <- model_theta_tau() +
    feat_selection(par_const(500), par_const(250),
                   population = 1, at_time = 1.7)

  set.seed(6688)
  sum_stats <- msmsSimFunc(dm, c(1, 5))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)

  set.seed(6688)
  sum_stats2 <- msmsSimFunc(dm, c(1, 5))
  expect_equal(sum_stats, sum_stats2)
})


test_that("msmsSimFunc works with inter-locus variation", {
  if (!checkForMsms(FALSE, TRUE)) skip('msms not installed')

  dm_tmp <- CoalModel(5, 2) +
    feat_mutation(par_range('theta', 1, 5), variance = 17) +
    sumstat_seg_sites()
  expect_true(hasInterLocusVariation(dm_tmp))

  set.seed(1100)
  sum_stats <- msmsSimFunc(dm_tmp, c(3))
  expect_is(sum_stats$seg_sites, 'list')
  expect_equal(length(sum_stats$seg_sites), 2)
})


test_that('simulating unphased data works', {
  model <- model_theta_tau() + feat_unphased(2, 1) + sumstat_seg_sites()
  stats <- msmsSimFunc(model, c(1,5))
  expect_equal(dim(stats$jsfs), c(11, 16))
  expect_equal(nrow(stats$seg_sites[[1]]), 25)

  model <- model_theta_tau() + feat_unphased(3, 2) + sumstat_seg_sites()
  stats <- msmsSimFunc(model, c(1,5))
  expect_equal(dim(stats$jsfs), c(21, 31))
  expect_equal(nrow(stats$seg_sites[[1]]), 50)
})
