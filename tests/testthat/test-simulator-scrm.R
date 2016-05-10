context("Simulator scrm")

test_that("printing command works", {
  scrm <- get_simulator("scrm")
  model <- model_theta_tau()
  expect_equal(scrm$get_cmd(model),
               "scrm 25 10 -I 2 10 15 -ej tau 2 1 -t theta ")
})
