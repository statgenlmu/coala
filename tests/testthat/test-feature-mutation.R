context('Feature Mutation')


test_that('Creation of simple mutation features works', {
  feat <- feat_mutation(par_const(5))
  expect_equal(get_feature_table(feat)$parameter[1], '5')
})


test_that('Creation of finite sites models works', {
  for (model in sg_mutation_models) {
    feat <- feat_mutation(par_const(5), model=model)
    expect_equal(get_feature_table(feat)$parameter, c('5', model))
  }
  expect_error( feat_mutation(par_const(5), model="BLUB"))

  # tstv-ratio
  feat <- feat_mutation(par_const(5), model='HKY', tstv_ratio = 2)
  expect_equal(get_feature_table(feat)$parameter, c('5', 'HKY', '2'))
  feat <- feat_mutation(par_const(5), model='F84', tstv_ratio = 2)
  expect_equal(get_feature_table(feat)$parameter, c('5', 'F84', '2'))
  expect_error( feat_mutation(par_const(5), model='GTR', tstv_ratio = 2))
  expect_error( feat_mutation(par_const(5), tstv_ratio = 2))

  # base frequencies
  feat <- feat_mutation(par_const(5), model='HKY', tstv_ratio = 2,
                        base_frequencies = c(.1, .2, .3, .4))
  expect_equal(get_feature_table(feat)$parameter,
               c('5', 'HKY', '2', '0.1', '0.2', '0.3', '0.4'))
  feat <- feat_mutation(par_const(5), model='F84', tstv_ratio = 2,
                        base_frequencies = c(.1, .2, .3, .4))
  expect_error(feat_mutation(par_const(5), model='GTR',
                             base_frequencies=c(.1, .2, .3, .4)))
  expect_error(feat_mutation(par_const(5), model='HKY',
                             base_frequencies=c(.1, .2, .3, .6)))

  # GTR rates
  feat <- feat_mutation(par_const(5), model='GTR', gtr_rates=1:6)
  expect_equal(get_feature_table(feat)$parameter, c('5', 'GTR', 1:6))
  expect_error(feat_mutation(par_const(5), model='GTR', gtr_rates=1:5))
})
