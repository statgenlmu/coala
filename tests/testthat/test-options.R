context('Package Options')

test_that('get and set msms\' path works', {
  set_msms_path('blub/msms')
  expect_equal(get_msms_path(), 'blub/msms')
})


test_that('get and set seqgen\' path works', {
  set_msms_path('blub/seqgen')
  expect_equal(get_msms_path(), 'blub/seqgen')
})
