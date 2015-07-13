context("Feature Selection")

test_that("generation of selection cmd works", {
  if (!has_msms()) skip("msms not installed")
  msms <- get_simulator("msms")
  model  <- model_theta_tau() +
    feat_selection(111, 222, population = 1, time = 5)
  cmd <- msms$get_cmd(model)
  expect_true(grepl("-N 10000", cmd))
  expect_true(grepl("-SI 5 2 5e-04 0", cmd))
  expect_true(grepl("-SAA 111", cmd))
  expect_true(grepl("-SAa 222", cmd))
  expect_true(grepl(" $", cmd))

  model <- model_theta_tau() +
    feat_selection(strength_A = 123, population = 1, time = 5)
  cmd <- msms$get_cmd(model)
  expect_true(grepl("-N 10000", cmd))
  expect_true(grepl("-SI 5 2 5e-04 0", cmd))
  expect_true(grepl("-SA 123", cmd))
  expect_true(grepl(" $", cmd))
})


test_that("msms can simulate selection", {
  if (!has_msms()) skip("msms not installed")
  # With one population
  model <- coal_model(5, 1, 100) + par_named("s") + par_named("t") +
    feat_selection(par_expr(2 * s), par_expr(.5 * s), time = par_expr(t * 2)) +
    feat_mutation(5) +
    sumstat_sfs()

  expect_equal(select_simprog(model)$get_name(), "msms")
  stat <- simulate(model, pars = c(s = 20, t = 0.01))
  expect_that(stat$sfs, is_a("numeric"))

  # With three population
  if (!has_msms()) skip("msms not installed")
  model <- coal_model(c(5, 5, 5), 1, 100) +
    feat_selection(1000, 500, 1, time = 0.01) +
    feat_mutation(5) +
    feat_migration(1, symmetric = TRUE) +
    sumstat_sfs()
  expect_equal(select_simprog(model)$get_name(), "msms")
  stat <- simulate(model)
  expect_that(stat$sfs, is_a("numeric"))

  # With additive selection
  model <- coal_model(5, 1, 100) +
    feat_selection(strength_A = 123, time = 0.03) +
    feat_mutation(5) +
    sumstat_sfs()
  expect_equal(select_simprog(model)$get_name(), "msms")
  stat <- simulate(model)
  expect_that(stat$sfs, is_a("numeric"))
})
