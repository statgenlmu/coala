context('Feature Growth')


test_that("generating ms cmd for growth works", {
  ms <- get_simulator("ms")
  model <- coal_model(4, 1) + feat_mutation(1)
  expect_equal(ms$get_cmd(model + feat_growth(5, 1)), "ms 4 1 -t 1 -G 5 ")
  expect_equal(ms$get_cmd(model + feat_growth(5, 1, 2)), "ms 4 1 -t 1 -eG 2 5 ")

  model <- coal_model(4:5, 1) + feat_mutation(1) +
    feat_migration(1, symmetric = TRUE)

  expect_equal(ms$get_cmd(model + feat_growth(5, 1)),
               "ms 9 1 -I 2 4 5 -t 1 -eM 0 1 -g 1 5 ")
  expect_equal(ms$get_cmd(model + feat_growth(5, 1, 1.5)),
               "ms 9 1 -I 2 4 5 -t 1 -eM 0 1 -eg 1.5 1 5 ")

  model <- model +
    par_range("a", 0, 1) +
    feat_growth(par_expr(2 * a), 1, par_expr(3 * a))
  expect_equal(ms$get_cmd(model),
               "ms 9 1 -I 2 4 5 -t 1 -eM 0 1 -eg 3 * a 1 2 * a ")
  expect_that(ms$simulate(model, .5), is_a('list'))
})
