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

