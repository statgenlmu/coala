context("Simulation Program seqgen")

test_that("test.generateSeqgenOptions", {
  if (!checkForSeqgen(FALSE, TRUE)) skip('seqgen not installed')
  dm.hky <- model_hky()
  opts <- generateSeqgenOptions(dm.hky, c(1, 10), 1, c(0, 0, 10, 0, 0), 1)
  opts <- strsplit(opts, " ")[[1]]
  expect_true("-l" %in% opts)
  expect_true("-p" %in% opts)
  expect_true("-z" %in% opts)
  expect_true("-q" %in% opts)
  expect_true("-mHKY" %in% opts)
  expect_true("-t" %in% opts)
  expect_true("-s" %in% opts)
})


test_that("generation of tree models works", {
  if (!checkForSeqgen(FALSE, TRUE)) skip('seqgen not installed')
  for (dm in list(model_hky(), model_f84(), model_gtr())) {
    dm.ms <-generateTreeModel(dm, 1)
    sum.stats <- simulate(dm.ms, pars=c(1, 5))
    expect_false(is.null(sum.stats$file))
    expect_true(file.exists(sum.stats$file[[1]]))
    unlink(sum.stats$file)
  }
})


test_that("test.seqgenSingleSimFunc", {
  if (!checkForSeqgen(FALSE, TRUE)) skip('seqgen not installed')
  seqgenSingleSimFunc = getSimProgram("seq-gen")$sim_func

  set.seed(100)
  sum.stats <- seqgenSingleSimFunc(model_hky(), c(1, 10))
  expect_true(is.list(sum.stats))
  expect_true(is.array(sum.stats$jsfs))
  expect_true(sum(sum.stats$jsfs) > 0)

  set.seed(100)
  sum.stats2 <- seqgenSingleSimFunc(model_hky(), c(1, 10))
  expect_equal(sum.stats2$jsfs, sum.stats$jsfs)
})


test_that("All example models can be simulated", {
  if (!checkForSeqgen(FALSE, TRUE)) skip('seqgen not installed')
  set.seed(12)
  for (model in list(model_hky(), model_f84(), model_gtr())) {
    sum_stats <- simulate(model, pars=c(1, 5))
    expect_true(sum(sum_stats$jsfs) > 0)
  }
})


test_that("test.RateHeterogenity", {
  skip("Temporarily deactivated")
  if (!checkForSeqgen(FALSE, TRUE)) skip('seqgen not installed')
  set.seed(12)
  dm.rh <- dm.addMutationRateHeterogenity(dm.hky, 0.1, 5, categories.number = 5)
  jsfs <- simulate(dm.rh, c(1, 10, 1))
  expect_true(sum(jsfs$jsfs) > 0)
})


test_that("test.seqgenWithMsms", {
  if (!checkForSeqgen(FALSE, TRUE)) skip('seqgen not installed')
  if (!checkForMsms(FALSE, TRUE)) skip('msms not installed')

  m1 <- model_hky() + feat_selection(500, 250, 1, 0.1)
  set.seed(4444)
  sum.stats <- simulate(m1, pars = c(1, 5))
  expect_false(is.null(sum.stats$jsfs))

  set.seed(4444)
  sum.stats2 <- simulate(m1, pars = c(1, 5))
  expect_equal(sum.stats2, sum.stats)
})


test_that("seq-gen can simulate trios", {
  if (!checkForSeqgen(FALSE, TRUE)) skip('seqgen not installed')
  dm.lt <- model_f84() +
    locus_trio(locus_length = c(10, 20, 10), distance = c(5, 5), group = 1) +
    locus_trio(locus_length = c(20, 10, 15), distance = c(7, 5), group = 1) +
    sumstat_seg_sites()

  sum.stats <- simulate(dm.lt, pars=c(1, 10))
  expect_that(sum(sum.stats$jsfs), is_less_than(sum(sapply(sum.stats$seg.sites, ncol))))
})


test_that("Error is thrown without an outgroup", {
  if (!checkForSeqgen(FALSE, TRUE)) skip('seqgen not installed')
  model <- CoalModel(c(3, 3), 10) +
    feat_mutation(par_range('theta', 5, 10), model = 'HKY') +
    feat_pop_merge(par_range('tau', .5, 1), 2, 1) +
    sumstat_jsfs()
  expect_error(simulate(model, pars = c(7.5, .75)))
})


test_that("Simulation of trios with unequal mutation rates works", {
  skip("Temporarily deactivated")
#   if (!test_seqgen) skip('seq-gen not installed')
#   if (!isJaathaVariable('seqgen.exe')) setJaathaVariable('seqgen.exe', 'seq-gen')
#   dm <- dm.setTrioMutationRates(dm_trios, '17', 'theta', group=2)
#   expect_equal(getThetaName(dm, outer = FALSE, group = 2), '17')
#   expect_equal(getThetaName(dm, outer = TRUE, group = 2), 'theta')
#   expect_equal(getThetaName(dm, outer = FALSE, group = 1), 'theta')
#
#   expect_equal(getThetaName(dm_trios, outer = FALSE, group = 2), 'theta')
#   expect_equal(getThetaName(dm_trios, outer = TRUE, group = 2), 'theta')
#
#
#   # generateSeqgenOptionsCmd
#   dm <- dm.setTrioMutationRates(dm_trios, 'BLUB', 'BLA', group=2)
#   grp_mdl <- generateGroupModel(dm, 2)
#   cmd <- lapply(generateSeqgenOptionsCmd(grp_mdl), paste0, collapse = "")
#   expect_equal(length(cmd), 3)
#   expect_equal(grep("BLUB", cmd[[2]]), 1)
#   expect_equal(grep("theta", cmd[[2]]), integer(0))
#   expect_equal(grep("BLA", cmd[[1]]), 1)
#   expect_equal(grep("BLA", cmd[[3]]), 1)
#   expect_equal(grep("theta", cmd[[1]]), integer(0))
#   expect_equal(grep("theta", cmd[[3]]), integer(0))
#
#   # generateSeqgenOptions
#   dm <- dm.setTrioMutationRates(dm_trios, '1', '2', group=2)
#   grp_mdl <- generateGroupModel(dm, 2)
#   ll <- get_locus_length_matrix(grp_mdl)[1,]
#   cmds <- generateSeqgenOptions(grp_mdl, c(1, 5), locus = 1,
#                                 locus_lengths = ll, 1)
#   expect_equal(grep(as.character(2/ll[1]), cmds[[1]]), 1)
#   expect_equal(grep(as.character(1/ll[3]), cmds[[2]]), 1)
#   expect_equal(grep(as.character(2/ll[5]), cmds[[3]]), 1)
#
#   ll <- get_locus_length_matrix(grp_mdl)[2,]
#   cmds <- generateSeqgenOptions(grp_mdl, c(1, 5), locus = 2,
#                                 locus_lengths = ll, 1)
#   expect_equal(grep(as.character(2/ll[1]), cmds[[1]]), 1)
#   expect_equal(grep(as.character(1/ll[3]), cmds[[2]]), 1)
#   expect_equal(grep(as.character(2/ll[5]), cmds[[3]]), 1)
#
#   # simulate
#   grp_mdl <- grp_mdl + sumstat_seg_sites()
#   ss <- simulate(grp_mdl, c(1,5))
#   expect_false(is.null(ss$seg.sites))
})


test_that('a more complicated model works', {
  model <- CoalModel(c(5,5,2), 1, 100) +
    feat_mutation(par_range('theta', .1, 40), model = 'HKY',
                  base_frequencies = c(0.26, 0.20, 0.22, 0.32),
                  tstv_ratio = 1.26) +
    feat_migration(par_range('m12', 0.001, 5), 1, 2) +
    feat_migration(par_range('m21', 0.001, 5), 2, 1) +
    feat_size_change(par_range('q', 0.05, 40), population = 2, at.time = 0) +
    par_range("s1", 0.01, 2) + par_range("s2", 0.01, 2) +
    feat_growth(par_expr(log(1/s1)/tau), population = 1, at.time = 0) +
    feat_growth(par_expr(log(q/s2)/tau), population = 2, at.time = 0) +
    feat_size_change(par_expr(s1+s2), population = 1, at.time = par_expr(tau)) +
    feat_pop_merge(par_range('tau', 0.001, 5), 2, 1) +
    feat_pop_merge(par_expr(2*tau), 3, 1) +
    feat_recombination(par_const(10)) +
    feat_outgroup(3) +
    sumstat_jsfs()

  stat <- simulate(model, pars = c(10, 0.5, 0.6, 4.5, 0.2, 0.1, 0.4))
  expect_true(sum(stat$jsfs) > 0)
})
