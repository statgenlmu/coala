context('Simulator Class')

test_that('Simulator base class works', {
  sim <- Simulator$new()
  expect_equal(sim$get_name(), 'TEMPLATE')
  expect_equal(sim$get_features(), NULL)
  expect_equal(sim$get_sumstats(), NULL)
  expect_error(sim$simulate())
  expect_error(sim$get_cmd())
  expect_true(is_simulator(sim))
})


test_that('simulators can be registered', {
  register_simulator(Simulator)
  expect_true(is_simulator(get_simulator("TEMPLATE")))

  expect_error(register_simulator(1:10))
  expect_error(register_simulator(Feature))
})


test_that("reducing of simulation commands works", {
  df <- data.frame(locus_number = rep(1, 10),
                   command = as.character(1:10))
  expect_equal(df, reduce_sim_commands(df))

  df$command <- c(rep("a", 5), rep("b", 5))
  expect_equivalent(reduce_sim_commands(df),
                    data.frame(locus_number = c(5, 5),
                               command = c("a", "b"),
                               stringsAsFactors = FALSE))
})
