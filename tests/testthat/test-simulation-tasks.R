context("Simulation Tasks")

test_that("simulation tasks can be created", {
  task <- create_sim_task(get_simulator("scrm"), 10, test_opt = 1:5)
  expect_equal(task$locus_number, 10)
  expect_equal(task$get_arg("test_opt"), 1:5)
  expect_error(task$get_arg("not_set"))
  expect_identical(task$get_simulator(), get_simulator("scrm"))
})


test_that("hashes of simulation tasks ignore the locus number", {
  task1 <- create_sim_task(get_simulator("scrm"), 10, test_opt = 1:5)
  task2 <- create_sim_task(get_simulator("scrm"), 20, test_opt = 1:5)
  expect_is(task1$hash(), "character")
  expect_equal(task1$hash(), task2$hash())
})


test_that("hashes of simulation tasks respect options", {
  task1 <- create_sim_task(get_simulator("scrm"), 10, a = 1:5)
  task2 <- create_sim_task(get_simulator("scrm"), 10, a = 1:4)
  task3 <- create_sim_task(get_simulator("scrm"), 10, a = 1:5, b = 2)
  expect_true(task1$hash() != task2$hash())
  expect_true(task1$hash() != task3$hash())
})


test_that("hashes of simulation tasks respect the simulator", {
  task1 <- create_sim_task(get_simulator("scrm"), 10, a = 1)
  task2 <- create_sim_task(test_simulator(), 10, a = 1)
  expect_true(task1$hash() != task2$hash())
})


test_that("reduce_sim_tasks reduces simulation tasks", {
  tasks <- list(create_sim_task(test_simulator(), 10, a = 1),
                create_sim_task(test_simulator(), 10, b = 1))
  expect_equal(reduce_sim_tasks(tasks), tasks)

  tasks <- list(create_sim_task(test_simulator(), 10, a = 1),
                create_sim_task(test_simulator(), 15, b = 1),
                create_sim_task(test_simulator(), 10, a = 1))
  tasks_reduced <- reduce_sim_tasks(tasks)
  expect_equal(length(tasks_reduced), 2)
  expect_equal(tasks_reduced[[1]]$locus_number, 20)
  expect_equal(tasks_reduced[[1]]$hash(), tasks[[1]]$hash())
  expect_equal(tasks_reduced[[2]]$locus_number, 15)
  expect_equal(tasks_reduced[[2]]$hash(), tasks[[2]]$hash())
})


test_that("generate_sim_tasks generates simple simulation tasks", {
  scrm <- get_simulator("scrm")
  model <- model_theta_tau()
  pars <- c(tau = 1, theta = 2)

  expect_equal(generate_sim_tasks(model, pars),
               list(scrm$create_task(model, pars, 10, 1)))

  model <- model + locus_single(1)
  expect_equal(generate_sim_tasks(model, pars),
               list(scrm$create_task(create_group_model(model, 1), pars, 10),
                    scrm$create_task(create_group_model(model, 2), pars, 1)))
})


test_that("combine_results combines results from different tasks", {
  model <- model_theta_tau() + sumstat_trees()
  scrm <- get_simulator("scrm")
  task <- scrm$create_task(model, c(tau = 1, theta = 5), 2)
  res1 <- scrm$simulate(model, task)
  res2 <- scrm$simulate(model, task)
  res3 <- scrm$simulate(model, task)
  results <- combine_results(list(res1, res2, res3))

  expect_equal(results$seg_sites, list(res1$seg_sites[[1]],
                                       res1$seg_sites[[2]],
                                       res2$seg_sites[[1]],
                                       res2$seg_sites[[2]],
                                       res3$seg_sites[[1]],
                                       res3$seg_sites[[2]]))
  expect_equal(results$trees, list(res1$trees[[1]],
                                   res1$trees[[2]],
                                   res2$trees[[1]],
                                   res2$trees[[2]],
                                   res3$trees[[1]],
                                   res3$trees[[2]]))
  expect_equal(results$cmds, list(res1$cmd, res2$cmd, res3$cmd))
  expect_identical(results$simulators, list(res1$simulator,
                                            res2$simulator,
                                            res3$simulator))
})
