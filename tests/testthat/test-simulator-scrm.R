context("Simulator scrm")

model <- coal_model(10, 10) +
  feat_mutation(5) +
  feat_recombination(10)
scrm <- get_simulator("scrm")


test_that("scrm cmds are generated", {
  cmd_tmpl <- scrm$create_cmd_template(model)
  expect_true(grepl("-t", cmd_tmpl))
  expect_true(grepl("-r", cmd_tmpl))
})


test_that("simulation with scrm work", {
  sim_task <- generate_sim_tasks(model, NULL)
  res <- scrm$simulate(model, sim_task[[1]])
  expect_gt(sum(res$seg_sites[[1]]), 0)
  expect_equal(length(res$seg_sites), 10)
})
