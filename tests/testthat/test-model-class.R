context("Model Class")


test_that("creating models works", {
  model <- coal_model(11:12, 111, 1234)
  expect_true(is.model(model))
  expect_equal(get_sample_size(model), 11:12)
  expect_equal(get_locus_number(model), 111)
  expect_equal(get_locus_length(model), 1234)

  expect_is(model$id, 'character')
})


test_that("adding parameters works", {
  model <- coal_model(11:12, 100) + par_range("p1", 1, 5)
  par_table <- get_parameter_table(model)
  expect_equal("p1", par_table$name)
  expect_equal(1, par_table$lower.range)
  expect_equal(5, par_table$upper.range)

  expect_that(get_parameter(model), is_a("list"))
  expect_equal(length(get_parameter(model)), 1)
  expect_true(is.par(get_parameter(model)[[1]]))

  model <- model + par_range("p2", 1, 5)
  expect_equal(length(get_parameter(model)), 2)

  test <- list (1:10)
  class(test) <- 'CSR_OBJ'
  expect_error(model + test)
})


test_that("adding features works", {
  dm <- coal_model(11:12, 100)
  dm <- dm + Feature$new('blub', 5)
  expect_equal(nrow(search_feature(dm, 'blub')), 1)
  expect_that(length(get_features(dm)), is_more_than(0))

  dm <- coal_model(11:12, 100)
  dm <- dm + Feature$new('bli', par_range('bla', 1, 5))
  expect_equal(nrow(search_feature(dm, 'bli')), 1)
  expect_true('bla' %in% get_parameter_table(dm)$name)

  dm <- coal_model(11:12, 100)
  dm <- dm + Feature$new('bli', par_range('bla', 1, 5), variance='15')
  expect_true(has_inter_locus_var(dm))

  dm <- coal_model(11:12, 100)
  dm <- dm + Feature$new('bli', par_range('bla', 1, 5), group=1, variance='15')
  expect_false(has_inter_locus_var(dm, 0))
  expect_true(has_inter_locus_var(dm, 1))
})


test_that("test get_summary_statistics", {
  expect_equal(get_summary_statistics(coal_model(1:2, 2)), character())
  expect_equal(get_summary_statistics(model_theta_tau()),  "jsfs")

  expect_equal(get_summary_statistics(model_grps(), 1), "jsfs")
  expect_equal(get_summary_statistics(model_grps(), 2), "jsfs")
  expect_equal(get_summary_statistics(model_grps(), 3), "jsfs")
})


test_that("generation of group models", {
  dm <- get_group_model(model_theta_tau(), 1)
  expect_equal(get_feature_table(dm), get_feature_table(model_theta_tau()))
  expect_equal(dm$sum_stats, model_theta_tau()$sum_stats)

  dm <- model_theta_tau() + locus_averaged(10, 23, group = 1)
  dm <- get_group_model(dm, 1)
  expect_equal(nrow(get_feature_table(dm)),
               nrow(get_feature_table(model_theta_tau())))
  expect_true(all(get_feature_table(dm)$group == 0))
  expect_equal(get_locus_length(dm), 23)

  dm.3 <- model_theta_tau() +
    locus_averaged(5, 23, group = 1) +
    locus_averaged(31, 30, group = 2)
  dm <- get_group_model(dm.3, 1)
  expect_equal(nrow(get_feature_table(dm)),
               nrow(get_feature_table(model_theta_tau())))
  expect_true(all(get_feature_table(dm)$group == 0))
  expect_equal(get_locus_length(dm), 23)
  dm <- get_group_model(dm.3, 2)
  expect_equal(nrow(get_feature_table(dm)),
               nrow(get_feature_table(model_theta_tau())))
  expect_true(all(get_feature_table(dm)$group == 0))
  expect_equal(get_locus_length(dm), 30)
  expect_equal(get_locus_number(dm), 31)

  dm <- model_theta_tau() + sumstat_seg_sites(group = 2)
  dm_1 <- get_group_model(dm, 1)
  expect_equal(get_summary_statistics(dm, 1), get_summary_statistics(dm_1))
  dm_2 <- get_group_model(dm, 2)
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


test_that("getting the Theta Name works", {
  expect_equal(get_mutation_par(model_theta_tau()), "theta")
  expect_equal(get_mutation_par(model_hky()), "theta")
  expect_equal(get_mutation_par(model_f84()), "theta")

  dm.test <- coal_model(11:12, 100) +
    feat_mutation(par_range('abcd', 1, 5))
  expect_equal(get_mutation_par(dm.test), "abcd")
})


test_that("test.parInRange", {
  check_par_range(model_theta_tau(), c(1, 5))
  check_par_range(model_theta_tau(), c(2, 7))
  check_par_range(model_theta_tau(), c(0.5, 7.7))
  expect_error(check_par_range(model_theta_tau(), c(0, 5)))
  expect_error(check_par_range(model_theta_tau(), c(0, -1)))
  expect_error(check_par_range(model_theta_tau(), c(10, 1)))
  expect_error(check_par_range(model_theta_tau(), c(100, 100)))
  expect_error(check_par_range(model_theta_tau(), 1))
  expect_error(check_par_range(model_theta_tau(), matrix(1, 2, 2)))
  expect_error(check_par_range(model_theta_tau(), NULL))
})


test_that("test that scaling of model works", {
  model <- coal_model(11:12, 10) +
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


test_that("searching features works", {
  model <- model_theta_tau()
  expect_equal(nrow(search_feature(model, type = "sample")), 2)
  expect_equal(nrow(search_feature(model, type = c("pop_merge", "sample"))), 3)
  expect_equal(nrow(search_feature(model, type = c("pop_merge", "sample"),
                                  pop.sink = NA)), 2)
  expect_equal(nrow(search_feature(model, time.point = "tau")), 1)
  expect_equal(nrow(search_feature(model, time.point = NA)), 1)

  expect_equal(nrow(search_feature(model, type = "sample")), 2)
})


test_that("get loci length and number works", {
  model <- coal_model(10, 11, 101) + locus_averaged(12, 102, group = 2)
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


test_that("Adding and Getting inter locus variation works", {
  expect_false(has_inter_locus_var(model_theta_tau()))

  dm_tmp <- model_theta_tau() + feat_recombination(5, variance = 3)
  expect_true(has_inter_locus_var(dm_tmp))

  dm_tmp <- model_theta_tau() + feat_recombination(5, group = 2, variance = 3)
  expect_false(has_inter_locus_var(dm_tmp))
  expect_false(has_inter_locus_var(dm_tmp, 1))
  expect_true(has_inter_locus_var(dm_tmp, 2))
})


test_that('setTrioMutationsRates works', {
  warning("test about model with trio mutation rates deactivated")
#   dm <- dm.setTrioMutationRates(dm_trios, '17', 'theta', group=2)
#   expect_equal(nrow(search_feature(dm, 'mutation', group=2)), 1)
#   expect_equal(search_feature(dm, 'mutation', group=2)$parameter, "17")
#   expect_equal(nrow(search_feature(dm, 'mutation_outer', group=2)), 1)
#   expect_equal(search_feature(dm, 'mutation_outer', group=2)$parameter,
#                "theta")
})


test_that('getting the available Populations works', {
  dm <- coal_model(10:11, 100)
  expect_equal(get_populations(dm), 1:2)
  expect_equal(get_populations(model_theta_tau()), 1:2)
  expect_equal(get_populations(model_hky()), 1:3)
})


test_that('Outgroup setting and getting works', {
  model <- coal_model(1:4 * 2, 100)

  for (i in 1:4) {
    expect_equal(get_outgroup(model + feat_outgroup(i)), i)
    expect_equal(get_outgroup_size(model + feat_outgroup(i)), 2 * i)
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
  model <- coal_model(5:6) +
    locus_trio(locus_length = c(10, 30, 50), distance = c(20, 40)) +
    locus_trio(locus_length = c(50, 30, 10), distance = c(40, 20)) +
    locus_averaged(2, 100, group = 2)

  expect_equal(conv_middle_to_trio_pos(.5, model),
               c(45 / 150, 105 / 150))

  expect_equal(conv_middle_to_trio_pos(15, model, relative_in = FALSE),
               c(45 / 150, 105 / 150))

  expect_equal(conv_middle_to_trio_pos(.5, model, relative_out = FALSE),
               c(45, 105))

  expect_equal(conv_middle_to_trio_pos(10, model,
                                       relative_out = FALSE,
                                       relative_in = FALSE), c(40, 100))


  expect_equal(conv_middle_to_trio_pos(10, model, group = 2,
                                       relative_out = FALSE,
                                       relative_in = FALSE),
               c(10, 10))
  expect_equal(conv_middle_to_trio_pos(.5, model, group = 2), c(.5, .5))

  ss <- matrix(0, 6, 5)
  attr(ss, 'positions') <- c(0.1, 0.5, 0.2, 0.6, 0.5, 1)
  attr(ss, 'locus') <- rep(c(-1, 0, 1), each = 2)
  expect_equal(get_snp_positions(list(ss, ss), model),
               list(c(1, 5, 36, 48, 125, 150) / 150,
                    c(5, 25, 96, 108, 145, 150) / 150))
  expect_equal(get_snp_positions(list(ss, ss), model, relative = FALSE),
               list(c(1, 5, 36, 48, 125, 150),
                    c(5, 25, 96, 108, 145, 150)))

  ss <- matrix(0, 6, 5)
  attr(ss, 'positions') <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1)
  attr(ss, 'locus') <- rep(0, 6)
  expect_equal(get_snp_positions(list(ss, ss), model, group=2, relative=TRUE),
               list(c(.1, .3, .5, .7, .9, 1), c(.1, .3, .5, .7, .9, 1)))
  expect_equal(get_snp_positions(list(ss, ss), model, group=2, relative=FALSE),
               list(c(10, 30, 50, 70, 90, 100), c(10, 30, 50, 70, 90, 100)))
})


test_that('getting the ploidy and individuals works', {
  model <- model_theta_tau()
  expect_equal(get_ploidy(model), 1L)
  expect_equal(get_samples_per_ind(model), 1L)
  sample_size <- get_sample_size(model)
  expect_false(is_unphased(model))

  model <- model_theta_tau() + feat_unphased(4, 2)
  expect_equal(get_ploidy(model), 4L)
  expect_equal(get_samples_per_ind(model), 2L)
  expect_equal(get_sample_size(model), sample_size * 2)
  expect_equal(get_sample_size(model, TRUE), sample_size * 4)
  expect_true(is_unphased(model))
})


test_that('print works on models', {
  # Printing an empty model works
  out <- capture.output(print(coal_model()))
  expect_that(length(out), is_more_than(0))

  # Printing parameters works
  out <- capture.output(print(coal_model(5) + par_range("abc", 1, 5)))
  expect_that(length(grep("abc", out)), is_more_than(0))
})


test_that('getting par names works', {
  expect_equal(get_par_names(coal_model()), character(0))

  model <- coal_model() + par_range("a", 1, 2) + par_range("b", 2, 3)
  expect_equal(get_par_names(model), c("a", "b"))
})
