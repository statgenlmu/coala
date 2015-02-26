context('Locus Class')

test_that('Creating loci works', {
  locus <- Locus$new(1013)
  expect_true(is.locus(locus))
  expect_equal(locus$get_length(), 1013)
  expect_equal(locus$get_number(), 1)

  locus <- Locus$new(1014, 10, 'blub', 2)
  expect_true(is.locus(locus))
  expect_equal(locus$get_length(), 1014)
  expect_equal(locus$get_number(), 10)
  expect_equal(locus$get_name(), 'blub')
  expect_equal(locus$get_group(), 2)

  expect_error(Locus$new("abc", 10, 'blub', 2))
  expect_error(Locus$new(-5, 10, 'blub', 2))
  expect_error(Locus$new(10, '10', 'blub', 2))
  expect_error(Locus$new(10, -3, 'blub', 2))
  expect_error(Locus$new(10, 5, 17, 2))
  expect_error(Locus$new(10, 5, 17, -1))
  expect_error(Locus$new(10, 5, 17, 'bla'))
})


test_that('Adding a locus to a model works', {
  model2 <- CoalModel(10) + locus_averaged(10, 1010)
  expect_equal(get_locus_number(model2), 10)
  expect_equal(get_locus_length(model2), 1010)

  model2 <- model2 +
    locus_single(1017, 2) +
    locus_single(1018, 2) +
    locus_single(1019, 2)
  expect_equal(get_locus_number(model2), 10)
  expect_equal(get_locus_length(model2), 1010)
  expect_equal(get_locus_number(model2, 2), 3)
  expect_equal(get_locus_length(model2, 2), 1:3 + 1016)
})


test_that('Adding a locus trio works', {
  m <- CoalModel(10) +
    locus_trio(locus_names = c(left='a', middle='b', right='c'),
               locus_length = c(10, 30, 50),
               distance = c(20, 40)) +
    locus_trio(locus_names = c('a', 'b', 'c'),
               locus_length = c(middle=30, left=10, right=50),
               distance = c(20, 40)) +
    locus_trio(locus_names = c(right='c', left='a', middle='b'),
               locus_length = c(10, 30, 50),
               distance = c(middle_right=40, left_middle=20))

  expect_equivalent(get_locus_length_matrix(m),
                    matrix(1:5 * 10, 3, 5, byrow = TRUE))
})
