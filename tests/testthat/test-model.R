context("Model Class")


test_that("creating models works", {
  model <- CoalModel(11:12, 111, 1234)
  expect_true(is.model(model))
  expect_equal(get_sample_size(model), 11:12)
  expect_equal(get_locus_number(model), 111)
  expect_equal(get_locus_length(model), 1234)

  expect_is(model$id, 'character')
})


test_that("adding parameters works", {
  model <- CoalModel(11:12, 100) + par_range("theta", 1, 5)
  par_table <- get_parameter_table(model)
  expect_equal("theta", par_table$name)
  expect_equal(1, par_table$lower.range)
  expect_equal(5, par_table$upper.range)

  test <- list (1:10)
  class(test) <- 'CSR_OBJ'
  expect_error(model + test)
})


test_that("adding features works", {
  dm <- CoalModel(11:12, 100)
  dm <- dm + Feature$new('blub', 5)
  expect_equal(nrow(searchFeature(dm, 'blub')), 1)

  dm <- CoalModel(11:12, 100)
  dm <- dm + Feature$new('bli', par_range('bla', 1, 5))
  expect_equal(nrow(searchFeature(dm, 'bli')), 1)
  expect_true('bla' %in% get_parameter_table(dm)$name)

  dm <- CoalModel(11:12, 100)
  dm <- dm + Feature$new('bli', par_range('bla', 1, 5), variance='15')
  expect_true(hasInterLocusVariation(dm))

  dm <- CoalModel(11:12, 100)
  dm <- dm + Feature$new('bli', par_range('bla', 1, 5), group=1, variance='15')
  expect_false(hasInterLocusVariation(dm, 0))
  expect_true(hasInterLocusVariation(dm, 1))

  dm <- CoalModel(11:12, 100) + (Feature$new('blub', 5) + par_range('a', 1, 5))
  expect_equal(nrow(get_parameter_table(dm)), 1)

  dm <- CoalModel(11:12, 100) + (Feature$new('blub', 5) + SumStat$new('5'))
  expect_equal(get_summary_statistics(dm), '5')
})


test_that("test get_summary_statistics", {
  expect_equal(get_summary_statistics(CoalModel(1:2, 2)), character())
  expect_equal(get_summary_statistics(model_theta_tau()),  "jsfs")

  expect_equal(get_summary_statistics(model_grps(), 1), "jsfs")
  expect_equal(get_summary_statistics(model_grps(), 2), "jsfs")
  expect_equal(get_summary_statistics(model_grps(), 3), "jsfs")
})


test_that("generation of group models", {
  dm <- generateGroupModel(model_theta_tau(), 1)
  expect_equal(dm$features, model_theta_tau()$features)
  expect_equal(dm$sum_stats, model_theta_tau()$sum_stats)

  dm <- model_theta_tau() + locus_averaged(10, 23, group = 1)
  dm <- generateGroupModel(dm, 1)
  expect_equal(nrow(dm$features), nrow(model_theta_tau()$features))
  expect_true(all(dm$features$group == 0))
  expect_equal(get_locus_length(dm), 23)

  dm.3 <- model_theta_tau() +
    locus_averaged(5, 23, group = 1) +
    locus_averaged(31, 30, group = 2)
  dm <- generateGroupModel(dm.3, 1)
  expect_equal(nrow(dm$features), nrow(model_theta_tau()$features))
  expect_true(all(dm$features$group == 0))
  expect_equal(get_locus_length(dm), 23)
  dm <- generateGroupModel(dm.3, 2)
  expect_equal(nrow(dm$features), nrow(model_theta_tau()$features))
  expect_true(all(dm$features$group == 0))
  expect_equal(get_locus_length(dm), 30)
  expect_equal(get_locus_number(dm), 31)

  dm <- model_theta_tau() + sumstat_seg_sites(group = 2)
  dm_1 <- generateGroupModel(dm, 1)
  expect_equal(get_summary_statistics(dm, 1), get_summary_statistics(dm_1))
  dm_2 <- generateGroupModel(dm, 2)
  expect_equal(get_summary_statistics(dm, 2),
               get_summary_statistics(dm_2, 'all'))
})

test_that("test.getGroups", {
  expect_equal(get_groups(model_theta_tau()), 1)
  dm <- model_theta_tau() + locus_averaged(10, 23, group = 1)
  expect_equal(get_groups(dm), 1)
  dm <- model_theta_tau() + locus_averaged(10, 32, group = 2)
  expect_equal(get_groups(dm), 1:2)
  dm <- dm + locus_averaged(10, 32, group = 4)
  expect_equal(get_groups(dm), c(1:2, 4))
})


test_that("test.getSampleSize", {
  expect_equal(get_sample_size(model_theta_tau()), c(10L, 15L))
})


test_that("test.getThetaName", {
  expect_equal(getThetaName(model_theta_tau()), "theta")
  expect_equal(getThetaName(model_hky()), "theta")
  expect_equal(getThetaName(model_f84()), "theta")

  dm.test <- CoalModel(11:12, 100) +
    feat_mutation(par_range('abcd', 1, 5))
  expect_equal(getThetaName(dm.test), "abcd")
})


test_that("test.parInRange", {
  checkParInRange(model_theta_tau(), c(1, 5))
  checkParInRange(model_theta_tau(), c(2, 7))
  checkParInRange(model_theta_tau(), c(0.5, 7.7))
  expect_error(checkParInRange(model_theta_tau(), c(0, 5)))
  expect_error(checkParInRange(model_theta_tau(), c(0, -1)))
  expect_error(checkParInRange(model_theta_tau(), c(10, 1)))
  expect_error(checkParInRange(model_theta_tau(), c(100, 100)))
  expect_error(checkParInRange(model_theta_tau(), 1))
  expect_error(checkParInRange(model_theta_tau(), matrix(1, 2, 2)))
  expect_error(checkParInRange(model_theta_tau(), NULL))
})


test_that("test that scaling of model works", {
  model <- CoalModel(11:12, 10) +
    locus_averaged(25, 10, group = 1) +
    locus_averaged(25, 15, group = 2) +
    locus_single(101, group = 3) +
    locus_single(102, group = 3)

  model <- scale_model(model, 5)
  expect_equal(get_locus_number(model, 0), 2L)
  expect_equal(get_locus_number(model, 1), 5L)
  expect_equal(get_locus_number(model, 2), 5L)
  expect_equal(get_locus_number(model, 3), 2L)
})


test_that("searchFeature", {
  model <- model_theta_tau()
  expect_equal(nrow(searchFeature(model, type = "sample")), 2)
  expect_equal(nrow(searchFeature(model, type = c("pop_merge", "sample"))), 3)
  expect_equal(nrow(searchFeature(model, type = c("pop_merge", "sample"),
                                  pop.sink = NA)), 2)
  expect_equal(nrow(searchFeature(model, time.point = "tau")), 1)
  expect_equal(nrow(searchFeature(model, time.point = NA)), 1)

  expect_equal(nrow(searchFeature(model, type = "sample")), 2)
})


test_that("get loci length and number works", {
  model <- CoalModel(10, 11, 101) + locus_averaged(12, 102, group = 2)
  expect_equal(get_locus_number(model), 11)
  expect_equal(get_locus_number(model, 1), 11)
  expect_equal(get_locus_number(model, 2), 12)
  expect_equal(get_locus_length(model, 1), 101)
  expect_equal(get_locus_length(model), 101)
  expect_equal(get_locus_length(model, 2), 102)

  model <- model + locus_averaged(21, 201, group = 1)
  expect_equal(get_locus_number(model), 21)
  expect_equal(get_locus_number(model, 1), 21)
  expect_equal(get_locus_number(model, 2), 12)
  expect_equal(get_locus_length(model, 1), 201)
  expect_equal(get_locus_length(model), 201)
  expect_equal(get_locus_length(model, 2), 102)
})


test_that('locus length matrix generations works', {
  # Multiple loci with equal length
  dimnames <- list(NULL,  c('length_l', 'length_il', 'length_m',
                            'length_ir', 'length_r') )

  expect_equal(get_locus_length_matrix(model_theta_tau()),
               matrix(c(0, 0, 1000, 0, 0), 10, 5, TRUE, dimnames))

  # Multiple loci with differnt length
  dm <- model_theta_tau() +
    locus_single(21, group = 2) +
    locus_single(22, group = 2) +
    locus_single(23, group = 2)

  expect_equal(get_locus_length_matrix(dm, group = 2),
               matrix(c(0, 0, 21, 0, 0,
                        0, 0, 22, 0, 0,
                        0, 0, 23, 0, 0), 3, 5, TRUE, dimnames))
})


test_that("simulation works", {
  expect_error(simulate(model_theta_tau(), pars=1))
  expect_error(simulate(model_theta_tau(), pars=1:3))
  expect_error(simulate(model_theta_tau(), pars=c(2, 50)))
  sum.stats <- simulate(model_theta_tau(), pars=c(1, 5))
  expect_true(is.list(sum.stats))
  expect_false(is.null(sum.stats$jsfs))
  expect_true(sum(sum.stats$jsfs) > 0)

  sum.stats <- simulate(model_grps(), pars=c(1, 5))
  expect_false(is.null(sum.stats))
  expect_false(is.null(sum.stats$jsfs.1))
  expect_true(sum(sum.stats$jsfs.1) > 0)
  expect_false(is.null(sum.stats$jsfs.2))
  expect_true(sum(sum.stats$jsfs.2) > 0)
  expect_false(is.null(sum.stats$jsfs.3))
  expect_true(sum(sum.stats$jsfs.3) > 0)
})


test_that("Adding and Getting inter locus variation works", {
  skip('model test for inter-locus-variation deactivated')
  expect_false(dm.hasInterLocusVariation(model_theta_tau()))

  dm_tmp <- dm.addInterLocusVariation(model_theta_tau(), 0)
  expect_true(dm.hasInterLocusVariation(dm_tmp))
  dm_tmp <- dm.addInterLocusVariation(dm_tmp, 0)
  expect_true(dm.hasInterLocusVariation(dm_tmp))
  expect_equal(nrow(searchFeature(dm_tmp, 'inter_locus_variation')), 1)

  dm_tmp <- dm.addInterLocusVariation(model_theta_tau(), 2)
  expect_false(dm.hasInterLocusVariation(dm_tmp))
  expect_false(dm.hasInterLocusVariation(dm_tmp, 1))
  expect_true(dm.hasInterLocusVariation(dm_tmp, 2))
})


test_that("test.printEmptyDM", {
  tmp_file <- tempfile()
  sink(tmp_file)
  dm <- CoalModel(25:26, 100)
  print(dm)
  sink(NULL)
  unlink(tmp_file)
})


test_that("test.printGroupDM", {
  tmp_file <- tempfile()
  sink(tmp_file)
  print(model_grps())
  sink(NULL)
  unlink(tmp_file)
})


test_that('setTrioMutationsRates works', {
  warning("test about model with trio mutation rates deactivated")
#   dm <- dm.setTrioMutationRates(dm_trios, '17', 'theta', group=2)
#   expect_equal(nrow(searchFeature(dm, 'mutation', group=2)), 1)
#   expect_equal(searchFeature(dm, 'mutation', group=2)$parameter, "17")
#   expect_equal(nrow(searchFeature(dm, 'mutation_outer', group=2)), 1)
#   expect_equal(searchFeature(dm, 'mutation_outer', group=2)$parameter, "theta")
})


test_that('getting the available Populations works', {
  dm <- CoalModel(10:11, 100)
  expect_equal(get_populations(dm), 1:2)
  expect_equal(get_populations(model_theta_tau()), 1:2)
  expect_equal(get_populations(model_hky()), 1:3)
})


test_that('Outgroup setting and getting works', {
  model <- CoalModel(1:4*2, 100)

  for (i in 1:4) {
    expect_equal(get_outgroup(model + feat_outgroup(i)), i)
    expect_equal(get_outgroup_size(model + feat_outgroup(i)), 2*i)
  }
})


test_that('get population individuals works', {
  expect_equal(get_population_indiviuals(model_theta_tau(), 1), 1:10)
  expect_equal(get_population_indiviuals(model_theta_tau(), 2), 11:25)
  expect_error(get_population_indiviuals(model_theta_tau(), 3))

  expect_equal(get_population_indiviuals(model_hky(), 1), 1:3)
  expect_equal(get_population_indiviuals(model_hky(), 2), 4:6)
  expect_equal(get_population_indiviuals(model_hky(), 3), 7)
})


test_that('converting positions for trios works', {
  model <- CoalModel(5:6) +
    locus_trio(locus_length = c(10, 30, 50), distance = c(20, 40)) +
    locus_trio(locus_length = c(50, 30, 10), distance = c(40, 20)) +
    locus_averaged(2, 100, group = 2)

  expect_equal(conv_middle_to_trio_pos(10, model, relative = FALSE),
               c(40, 100))
  expect_equal(conv_middle_to_trio_pos(.5, model, relative = TRUE),
               c(45/150, 105/150))

  expect_equal(conv_middle_to_trio_pos(10, model, group = 2, relative = FALSE),
               c(10, 10))
  expect_equal(conv_middle_to_trio_pos(.5, model, group = 2, relative = TRUE),
               c(.5, .5))


  ss <- matrix(0, 6, 5)
  attr(ss, 'positions') = c(0.1, 0.5, 0.2, 0.6, 0.5, 1)
  attr(ss, 'locus') = rep(c(-1, 0, 1), each = 2)
  expect_equal(get_snp_positions(list(ss, ss), model),
               list(c(1, 5, 36, 48, 125, 150) / 150,
                    c(5, 25, 96, 108, 145, 150) / 150))
  expect_equal(get_snp_positions(list(ss, ss), model, relative = FALSE),
               list(c(1, 5, 36, 48, 125, 150),
                    c(5, 25, 96, 108, 145, 150)))

  ss <- matrix(0, 6, 5)
  attr(ss, 'positions') = c(0.1, 0.3, 0.5, 0.7, 0.9, 1)
  attr(ss, 'locus') = rep(0, 6)
  expect_equal(get_snp_positions(list(ss, ss), model, group=2, relative=TRUE),
               list(c(.1, .3, .5, .7, .9, 1), c(.1, .3, .5, .7, .9, 1)))
  expect_equal(get_snp_positions(list(ss, ss), model, group=2, relative=FALSE),
               list(c(10, 30, 50, 70, 90, 100), c(10, 30, 50, 70, 90, 100)))
})
