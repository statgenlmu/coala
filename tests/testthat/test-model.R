context("Model Class")


test_that("creating models works", {
  model <- CoalModel(11:12, 111, 1234)
  expect_true(is.model(model))
  expect_equal(get_sample_size(model), 11:12)
  expect_equal(get_locus_number(model), 111)
  expect_equal(get_locus_length(model), 1234)
})


test_that("adding parameters works", {
  model <- CoalModel(11:12, 100) + par_range("theta", 1, 5)
  par_table <- get_parameter_table(model)
  expect_equal("theta", par_table$name)
  expect_equal(1, par_table$lower.range)
  expect_equal(5, par_table$upper.range)
  expect_error(model + 5)
  expect_error(model + 'bla')
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
})


test_that("test get_summary_statistics", {
  expect_equal(get_summary_statistics(CoalModel(1:2, 2)), character())
  expect_equal(get_summary_statistics(model_theta_tau()),  "jsfs")

  expect_equal(get_summary_statistics(model_grps(), 1), "jsfs")
  expect_equal(get_summary_statistics(model_grps(), 2), "jsfs")
  expect_equal(get_summary_statistics(model_grps(), 3), "jsfs")
})



test_that("test model finalization", {
  dm <- dm.finalize(model_grps())
  dm.1 <- dm$options$grp.models[["1"]]
  dm.2 <- dm$options$grp.models[["2"]]
  dm.3 <- dm$options$grp.models[["3"]]
  expect_false(is.null(dm.1))
  expect_false(is.null(dm.2))
  expect_false(is.null(dm.3))
  expect_false(is.null(dm.1$options$ms.cmd))
  expect_false(is.null(dm.2$options$ms.cmd))
  expect_false(is.null(dm.3$options$ms.cmd))
})


test_that("test.generateGroupModel", {
  dm <- generateGroupModel(model_theta_tau(), 1)
  expect_equal(dm$features, model_theta_tau()$features)
  expect_equal(dm$sum_stats, model_theta_tau()$sum_stats)
  expect_equal(dm$options, model_theta_tau()$options)

  dm <- dm.addLocus(model_theta_tau(), 23, 10, group = 1)
  dm <- generateGroupModel(dm, 1)
  expect_equal(nrow(dm$features), nrow(model_theta_tau()$features))
  expect_true(all(dm$features$group == 0))
  expect_equal(get_locus_length(dm), 23)

  dm.3 <- dm.addLocus(model_theta_tau(), 23, 5, group = 1)
  dm.3 <- dm.addLocus(dm.3, 30, 31, group = 2)
  dm <- generateGroupModel(dm.3, 1)
  expect_equal(nrow(dm$features), nrow(model_theta_tau()$features))
  expect_true(all(dm$features$group == 0))
  expect_equal(get_locus_length(dm), 23)
  dm <- generateGroupModel(dm.3, 2)
  expect_equal(nrow(dm$features), nrow(model_theta_tau()$features))
  expect_true(all(dm$features$group == 0))
  expect_equal(get_locus_length(dm), 30)
  expect_equal(get_locus_number(dm), 31)
  sum.stats <- model_theta_tau()$sum_stats

  dm <- model_theta_tau() + sumstat_seg_sites(group = 2)
  dm.1 <- generateGroupModel(dm, 1)
  expect_true(sum.stats$name %in% dm.1$sum_stats$name)
  expect_true(dm.1$sum_stats$name %in% sum.stats$name)
  dm.2 <- generateGroupModel(dm, 2)
  expect_equal(nrow(dm.2$sum_stats), nrow(sum.stats) + 1)
  expect_true("seg.sites" %in% get_summary_statistics(dm.2))
  dm$options[["bli"]] <- 0
  dm$options[["bla"]] <- 0
  dm$options[["blub"]] <- 0
  dm$options[["group.1"]] <- list(blub = 1, bli = 1)
  dm$options[["group.2"]] <- list(blub = 2, bla = 2)
  dm.1 <- generateGroupModel(dm, 1)
  dm.2 <- generateGroupModel(dm, 2)
  expect_equal(dm.1$options$bli, 1)
  expect_equal(dm.1$options$bla, 0)
  expect_equal(dm.1$options$blub, 1)
  expect_equal(dm.2$options$bli, 0)
  expect_equal(dm.2$options$bla, 2)
  expect_equal(dm.2$options$blub, 2)
})

test_that("test.getGroups", {
  expect_equal(get_groups(model_theta_tau()), 1)
  dm <- dm.addLocus(model_theta_tau(), 23, 10, group = 1)
  expect_equal(get_groups(dm), 1)
  dm <- dm.addLocus(model_theta_tau(), 32, 10, group = 2)
  expect_equal(get_groups(dm), 1:2)
  dm <- dm.addLocus(dm, 32, group = 4)
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


test_that("test.scaleDemographicModel", {
  dm <- CoalModel(11:12, 10)
  dm <- dm.addLocus(dm, 10, 25, group = 1)
  dm <- dm.addLocus(dm, 15, 25, group = 2)
  dm <- scaleDemographicModel(dm, 5)
  expect_equal(get_locus_number(dm, 0), 2L)
  expect_equal(get_locus_number(dm, 1), 5L)
  expect_equal(get_locus_number(dm, 2), 5L)
})


test_that("searchFeature", {
  expect_equal(nrow(searchFeature(model_theta_tau(), type = "sample")), 2)
  expect_equal(nrow(searchFeature(model_theta_tau(), type = c("pop_merge", "sample"))), 3)
  expect_equal(nrow(searchFeature(model_theta_tau(), type = c("pop_merge", "sample"),
                                  pop.sink = NA)), 2)
  expect_equal(nrow(searchFeature(model_theta_tau(), time.point = "tau")), 1)
  expect_equal(nrow(searchFeature(model_theta_tau(), time.point = NA)), 1)
})


test_that("set and get loci length works", {
  dm_tmp <- dm.setLociLength(model_theta_tau(), 17)
  expect_equal(get_locus_length(dm_tmp), 17)

  dm_tmp <- dm.addLocus(dm_tmp, 22, 10, group = 2)
  dm_tmp <- dm.setLociLength(dm_tmp, 23, group = 2)
  expect_equal(get_locus_length(dm_tmp), 17)
  expect_equal(get_locus_length(dm_tmp, 1), 17)
  expect_equal(get_locus_length(dm_tmp, 2), 23)

  dm_tmp <- dm.addLocus(dm_tmp, 32, 10, group = 1)
  dm_tmp <- dm.setLociLength(dm_tmp, 32, group = 1)
  expect_equal(get_locus_length(dm_tmp), 32)
  expect_equal(get_locus_length(dm_tmp, 1), 32)
  expect_equal(get_locus_length(dm_tmp, 2), 23)
})


test_that("set and get loci number works", {
  dm <- dm.setLociNumber(model_theta_tau(), 17)
  expect_equal(get_locus_number(dm), 17)

  dm <- dm.addLocus(dm, 22, 10, group = 2)
  dm <- dm.setLociNumber(dm, 23, group = 2)
  expect_equal(get_locus_number(dm), 17)
  expect_equal(get_locus_number(dm, 1), 17)
  expect_equal(get_locus_number(dm, 2), 23)

  dm <- dm.addLocus(dm, 32, 10, group = 1)
  dm <- dm.setLociNumber(dm, 32, group = 1)
  expect_equal(get_locus_number(dm), 32)
  expect_equal(get_locus_number(dm, 1), 32)
  expect_equal(get_locus_number(dm, 2), 23)
})


test_that('locus length matrix generations works', {
  # Multiple loci with equal length
  dimnames <- list(NULL,  c('length_l', 'length_il', 'length_m',
                            'length_ir', 'length_r') )

  expect_equal(get_locus_length_matrix(model_theta_tau()),
               matrix(c(0, 0, 1000, 0, 0), 10, 5, TRUE, dimnames))

  # Multiple loci with differnt length
  dm <- dm.addLocus(model_theta_tau(), 21, 1, group = 2)
  dm <- dm.addLocus(dm, 22, 1, group = 2)
  dm <- dm.addLocus(dm, 23, 1, group = 2)
  expect_equal(get_locus_length_matrix(dm, group = 2),
               matrix(c(0, 0, 21, 0, 0,
                        0, 0, 22, 0, 0,
                        0, 0, 23, 0, 0), 3, 5, TRUE, dimnames))
})


test_that("test simulation works", {
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


test_that("Loci trios are added to model", {
  dm.lt <- dm.addLocusTrio(model_hky(), locus_length =  c(3, 5, 7),
                           distance = c(4, 6), group = 2)

  expect_true(all(get_locus_length_matrix(dm.lt, 2) == 3:7))
  expect_true(nrow(searchFeature(dm.lt, 'locus_trios', group = 2)) > 0)
})


test_that("Adding and Getting inter locus variation works", {
  warning('model test for inter-locus-variation deactivated')
#   expect_false(dm.hasInterLocusVariation(model_theta_tau()))
#
#   dm_tmp <- dm.addInterLocusVariation(model_theta_tau(), 0)
#   expect_true(dm.hasInterLocusVariation(dm_tmp))
#   dm_tmp <- dm.addInterLocusVariation(dm_tmp, 0)
#   expect_true(dm.hasInterLocusVariation(dm_tmp))
#   expect_equal(nrow(searchFeature(dm_tmp, 'inter_locus_variation')), 1)
#
#   dm_tmp <- dm.addInterLocusVariation(model_theta_tau(), 2)
#   expect_false(dm.hasInterLocusVariation(dm_tmp))
#   expect_false(dm.hasInterLocusVariation(dm_tmp, 1))
#   expect_true(dm.hasInterLocusVariation(dm_tmp, 2))
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
