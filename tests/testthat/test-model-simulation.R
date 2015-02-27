context('Model Simulation')

test_that("models can be simulated", {
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


test_that("model without ranged parameters can be simualted", {
  model <- coal_model(5, 1) + feat_mutation(5) + sumstat_sfs("sfs")
  stat <- simulate(model)
  expect_false(is.null(stat$sfs))
  stat <- simulate(model, pars=c())
  expect_false(is.null(stat$sfs))
})
