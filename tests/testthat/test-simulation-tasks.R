context("Simulation Tasks")


test_that("generate_sim_tasks generates simple simulation tasks", {
  register_simulator(test_simulator())
  model <- test_model(locus_number = 10)

  expect_equal(generate_sim_tasks(model, parameters = c(a = 1, b = 2)),
               data.frame(simulator = "test",
                          locus_number = 10,
                          command = "sum( 1 + 2 )",
                          stringsAsFactors = FALSE))

  expect_equal(generate_sim_tasks(model + locus_single(1),
                                  parameters = c(a = 1, b = 2)),
               data.frame(simulator = c("test", "test"),
                          locus_number = c(10, 1),
                          command = rep("sum( 1 + 2 )", 2),
                          stringsAsFactors = FALSE))
})


test_that("combine_results combines results from different tasks", {
  model <- model_theta_tau() + sumstat_trees()
  scrm <- get_simulator("scrm")
  res1 <- scrm$simulate(model, 2, "-t 5 -T")
  res2 <- scrm$simulate(model, 1, "-r 1 10 -t 5 -T")
  res3 <- scrm$simulate(model, 1, "-r 2 10 -t 5 -T")
  results <- combine_results(list(res1, res2, res3))

  expect_equal(results$seg_sites, list(res1$seg_sites[[1]],
                                       res1$seg_sites[[2]],
                                       res2$seg_sites[[1]],
                                       res3$seg_sites[[1]]))
  expect_equal(results$trees, list(res1$trees[[1]],
                                   res1$trees[[2]],
                                   res2$trees[[1]],
                                   res3$trees[[1]]))
  expect_equal(results$cmds, list(res1$cmd, res2$cmd, res3$cmd))
  expect_identical(results$simulators, list(res1$simulator,
                                            res2$simulator,
                                            res3$simulator))
})


test_that("generate_sim_tasks generates special tasks for seqgen", {
  tasks <- generate_sim_tasks(model_hky(), parameters = c(tau = 1, theta = 5))
  tasks <- generate_sim_tasks(model_hky() + locus_trio(), parameters = c(tau = 1, theta = 5))
})
