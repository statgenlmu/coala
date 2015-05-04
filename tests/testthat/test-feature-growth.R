context('Feature Growth')

test_that('Creation of growth features works', {
  expect_equal(feat_growth(2, 1)$get_rate(), "2")
  expect_equal(feat_growth(2, 1)$get_population(), 1)

  feat <- feat_growth(3, 2, 5)
  expect_equal(feat$get_rate(), "3")
  expect_equal(feat$get_population(), 2)
  expect_equal(feat$get_time(), "5")

  expect_error(feat_growth(3, "A", 5))
  expect_error(feat_growth(3, 1:2, 5))
})


test_that("generating ms cmd for growth works", {
  model <- coal_model(15, 1) + feat_growth(par_range("alpha", 0, 1), 1)
  expect_equal(get_simulator("ms")$get_cmd(model), "ms 15 1 -eg 0 1 alpha ")
  model <- coal_model(15, 1) + feat_growth(5, 1)
  expect_equal(get_simulator("ms")$get_cmd(model), "ms 15 1 -eg 0 1 5 ")
})
