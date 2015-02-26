context("Feature Selection")

test_that("a variable start works", {
  model <- CoalModel(5) +
    feat_selection(strength_AA = par_const(0),
                   strength_Aa = par_const(1000),
                   population = 1,
                   at_time = par_range('t_sel', 0.001, 0.1))

  expect_equal(get_parameter_table(model), data.frame(name="t_sel",
                                                      lower.range=0.001,
                                                      upper.range=0.1,
                                                      stringsAsFactors = FALSE))
  expect_true(all(get_feature_table(model)$time.point %in% c('0', 't_sel')))

  model <- CoalModel(5) +
    feat_selection(strength_AA = par_const(0),
                   strength_Aa = par_const(1000),
                   population = 1,
                   at_time = par_expr(2 * x))

  expect_true(all(get_feature_table(model)$time.point %in% c('0', '2 * x')))

  model <- CoalModel(5) +
    feat_selection(strength_AA = par_const(0),
                   strength_Aa = par_const(1000),
                   population = 1,
                   at_time = 0.2)
  expect_true(all(get_feature_table(model)$time.point %in% c('0', '0.2')))
})
