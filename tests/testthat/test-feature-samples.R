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
  expect_equal(get_sample_size(coal_model(10)), 10)
})


test_that("generating ms cmd works", {
  expect_equal(conv_to_ms_arg(feat_sample(2), NULL), "")
  expect_equal(conv_to_ms_arg(feat_sample(1:2), NULL), "-I 2 1 2 ")
  expect_equal(conv_to_ms_arg(feat_sample(1:3), NULL), "-I 3 1 2 3 ")

  model <- coal_model(15, 1)
  expect_equal(get_simulator("ms")$get_cmd(model), "ms 15 1 ")
  model <- coal_model(15:16, 1)
  expect_equal(get_simulator("ms")$get_cmd(model), "ms 31 1 -I 2 15 16 ")
})
