context("Model Class")


test_that("creating models works", {
  model <- coal_model(11:12, 111, 1234)
  expect_true(is.model(model))
  expect_equal(get_sample_size(model), 11:12)
  expect_equal(get_locus_number(model), 111)
  expect_equal(get_locus_length(model, 1), 1234)

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
  class(test) <- 'BLUB'
  expect_error(model + test)
})


test_that("adding features works", {
  expect_equal(length(get_features(coal_model())), 0)
  expect_equal(nrow(get_feature_table(coal_model())), 0)

  dm <- coal_model(11:12, 100)
  dm <- dm + Feature$new('blub', 5)
  expect_equal(nrow(search_feature(dm, 'blub')), 1)
  expect_that(length(get_features(dm)), is_more_than(0))
  expect_true('blub' %in% get_feature_table(dm)$type)


  dm <- coal_model(11:12, 100)
  dm <- dm + Feature$new('bli', par_range('bla', 1, 5))
  expect_equal(nrow(search_feature(dm, 'bli')), 1)
  expect_true('bla' %in% get_parameter_table(dm)$name)

  dm <- coal_model(11:12, 100)
  dm <- dm + Feature$new('bli', par_range('bla', 1, 5), variance='15')
  expect_true(has_inter_locus_var(dm))
})


test_that("test get_summary_statistics", {
  expect_equal(get_summary_statistics(coal_model(1:2, 2)), character())
  expect_equal(get_summary_statistics(model_theta_tau()),  "jsfs")

  expect_equal(get_summary_statistics(model_grps(), 1), "jsfs")
  expect_equal(get_summary_statistics(model_grps(), 2), "jsfs")
  expect_equal(get_summary_statistics(model_grps(), 3), "jsfs")
})


test_that("sample sizes are reported corrently", {
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
  skip("temorarily deactivated")
  model <- coal_model(11:12, 10) +
    locus_averaged(25, 10) +
    locus_averaged(25, 15) +
    locus_single(101) +
    locus_single(102)

  model <- scale_model(model, 5)

  expect_equal(get_locus_number(model), 5L)
  expect_equal(get_locus_number(model), 5L)
  expect_equal(get_locus_number(model), 2L)
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
  expect_equal(nrow(search_feature(model, type = "mutation",
                                   pop.source = NA)), 1)
  expect_equal(nrow(search_feature(model, type = "pop_merge",
                                   pop.sink = 1)), 1)
})


test_that("get loci length and number works", {
  model <- coal_model(10, 11, 101) +
    locus_averaged(12, 102) +
    locus_trio(locus_length = 1:3, distance = 10:11)

  expect_equal(get_locus_number(model), 24)
  expect_equal(get_locus_length(model, 1), 101)
  expect_equal(get_locus_length(model, 5), 101)
  expect_equal(get_locus_length(model, 11), 101)
  expect_equal(get_locus_length(model, 15), 102)
  expect_equal(get_locus_length(model, 23), 102)
  expect_equal(get_locus_length(model, 24), 6)

  expect_equivalent(get_locus_length(model, 1, total = FALSE), 101)
  expect_equivalent(get_locus_length(model, 24, total = FALSE), 1:3)

  expect_error(get_locus_length(model))
})


test_that('locus length matrix generations works', {
  # Multiple loci with equal length
  expect_equivalent(get_locus_length_matrix(model_theta_tau()),
                    matrix(c(0, 0, 1000, 0, 0, 1), 10, 6, TRUE))
  expect_equivalent(get_locus_length_matrix(model_theta_tau(), FALSE),
                    matrix(c(0, 0, 1000, 0, 0, 10), 1, 6, TRUE))

  # Multiple loci with differnt length
  model <- model_theta_tau() +
    locus_single(21) +
    locus_single(22) +
    locus_single(23)

  expect_equivalent(get_locus_length_matrix(model, individual = FALSE),
                    matrix(c(0, 0, 1000, 0, 0, 10,
                             0, 0, 21, 0, 0, 1,
                             0, 0, 22, 0, 0, 1,
                             0, 0, 23, 0, 0, 1), 4, 6, TRUE))
})


test_that("Adding and Getting inter locus variation works", {
  expect_false(has_inter_locus_var(model_theta_tau()))

  dm_tmp <- model_theta_tau() + feat_recombination(5, variance = 3)
  expect_true(has_inter_locus_var(dm_tmp))
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

  # Printing loci works
  out <- capture.output(print(coal_model(5) + locus_single(3131)))
  expect_that(length(grep("3131", out)), is_more_than(0))
})


test_that('getting par names works', {
  expect_equal(get_par_names(coal_model()), character(0))

  model <- coal_model() + par_range("a", 1, 2) + par_range("b", 2, 3)
  expect_equal(get_par_names(model), c("a", "b"))
})
