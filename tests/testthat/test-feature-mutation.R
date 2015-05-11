context('Feature Mutation')


test_that('Creation of finite sites features works', {
  expect_error(feat_mutation(5, model = "BLUB"))
  expect_error(feat_mutation(5, model = "F84"))

  # HKY
  feat <- feat_mutation(17, 'HKY', tstv_ratio = 2,
                        base_frequencies = rep(.25, 4))
  expect_equal(feat$get_tstv_ratio(), 2)
  expect_equal(feat$get_base_frequencies(), rep(.25, 4))
  expect_error(feat_mutation(5, model = "HKY"))
  expect_error(feat_mutation(5, model = "HKY", tstv_ratio = 2))
  expect_error(feat_mutation(5, model = "HKY", base_frequencies = rep(.25, 4)))
  expect_error(feat_mutation(5, model = "HKY", tstv_ratio = 2,
                             base_frequencies = rep(.5, 4)))


  # GTR rates
  feat <- feat_mutation(5, model = 'GTR', gtr_rates = 1:6)
  expect_equal(feat$get_gtr_rates(), 1:6)
  expect_error(feat_mutation(5, model = 'GTR', gtr_rates = 1:5))
  expect_error(feat_mutation(5, model = 'GTR', gtr_rates = "Blub"))
})


test_that("Parsing mutation to ms args works", {
  feat <- feat_mutation(5)
  ms_arg <- conv_to_ms_arg(feat, NULL)
  expect_that(ms_arg, is_a("character"))
  expect_true(grepl("-t", ms_arg))
  expect_true(grepl("5", ms_arg))

  ms <- get_simulator("ms")
  model <- coal_model(15, 1) + feat_mutation(par_range("theta", 1, 2))
  expect_equal(ms$get_cmd(model), "ms 15 1 -t theta ")
  model <- coal_model(15, 1) + feat_mutation(5)
  expect_equal(ms$get_cmd(model), "ms 15 1 -t 5 ")
  model <- coal_model(15, 1) +
    par_range("x", 1, 2) +
    feat_mutation(par_expr(2 * x))
  expect_equal(ms$get_cmd(model), "ms 15 1 -t 2 * x ")

  expect_error(conv_to_ms_arg(feat_mutation(5, "GTR"), NULL))
})

test_that("Parsing mutation to seqgen args works", {
  sg <- get_simulator("seqgen")
  model <- coal_model(10, 1) + feat_mutation(5, model = 'GTR', gtr_rates = 1:6)
  sg$get_cmd(model)
})
