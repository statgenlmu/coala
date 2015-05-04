context('Feature Sample')


test_that('Creation of sample features works', {
  expect_equal(feat_sample(2)$get_sizes(), 2)
  expect_equal(feat_sample(1:5)$get_sizes(), 1:5)
  expect_error(feat_sample("blub"))
  expect_error(feat_sample(numeric(0)))
  expect_error(feat_sample(2, 5))
})


test_that("sample sizes are reported corrently", {
  expect_equal(get_sample_size(coal_model(c(10, 15))), c(10, 15))
})
