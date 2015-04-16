context('Model Simulation')

test_that("basic models can be simulated", {
  expect_error(simulate(model_theta_tau(), pars=1))
  expect_error(simulate(model_theta_tau(), pars=1:3))
  expect_error(simulate(model_theta_tau(), pars=c(2, 50)))
  sum.stats <- simulate(model_theta_tau(), pars=c(1, 5))
  expect_true(is.list(sum.stats))
  expect_false(is.null(sum.stats$jsfs))
  expect_true(sum(sum.stats$jsfs) > 0)
})


test_that("models with multiple loci groups can be simulated", {
  model <- model_theta_tau() + locus_averaged(2, 100) + sumstat_seg_sites()
  sum_stats <- simulate(model, pars=c(1, 5))
  expect_equal(length(sum_stats$seg_sites), get_locus_number(model))
})


test_that("model without ranged parameters can be simualted", {
  model <- coal_model(5, 1) + feat_mutation(5) + sumstat_sfs("sfs")
  stat <- simulate(model)
  expect_false(is.null(stat$sfs))
})
