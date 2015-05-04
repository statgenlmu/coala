context('Feature Outgroup')

test_that('Creation of outgroup feature works', {
  expect_equal(feat_outgroup(1)$get_population(), 1)
  expect_equal(feat_outgroup(2)$get_population(), 2)
})
