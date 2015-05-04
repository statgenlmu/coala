context('Feature Class')

test_that('Creating features works', {
  skip("Deativated")
  feat <- Feature$new('abc', 5)
  expect_true(is.feature(feat))
  expect_equal(feat$get_table(),
               create_feature_table('abc', '5', NA, NA, NA))
  expect_false(feat$get_inter_locus_var())
  expect_equal(feat$get_parameters(), list())

  feat <- Feature$new('a', par_expr(2 * x), 2, 1, 't5')
  expect_true(is.feature(feat))
  expect_equal(feat$get_table(),
               create_feature_table('a', '2 * x', 2, 1, 't5'))
  expect_false(feat$get_inter_locus_var())
  x <- 5
  expect_equal(feat$get_parameters()[["1"]]$eval(), 10)

  # Test with two parameters
  feat <- Feature$new('blub', parameter=par_range('a', 5, 7),
                      time_point=par_range('b', 1, 1.5))
  expect_equal(length(feat$get_parameters()), 2)
  expect_equal(feat$get_table(),
               create_feature_table('blub', 'a', NA, NA, 'b'))

  # Test variance
  feat <- Feature$new('mutation', par_range('theta', 5, 7), variance='10')
  expect_true(is.feature(feat))
  expect_true(feat$get_inter_locus_var())

  par_expr <- feat$get_table()$parameter
  theta <- 5
  expect_true(eval(parse(text=par_expr)) != 5)
  sim <- sapply(1:1000, function(x) eval(parse(text=par_expr)))
  expect_true(abs(mean(sim) - theta) < .3)

  # Test zero.inflation
  feat <- Feature$new('mutation', par_range('theta', 5, 7), zero_inflation='.1')
  expect_true(is.feature(feat))
  expect_true(feat$get_inter_locus_var())

  par_expr <- feat$get_table()$parameter
  locus_number <- 100
  locus <- 1; expect_equal(eval(parse(text=par_expr)), 0)
  locus <- 5; expect_equal(eval(parse(text=par_expr)), 0)
  locus <- 10; expect_equal(eval(parse(text=par_expr)), 0)
  locus <- 11; expect_equal(eval(parse(text=par_expr)), 5)
  locus <- 30; expect_equal(eval(parse(text=par_expr)), 5)
  locus <- 72; expect_equal(eval(parse(text=par_expr)), 5)

  # Test zero.inflation & variance
  set.seed(111222333)
  feat <- Feature$new('mutation', parameter = par_range('theta', 5, 7),
                      variance = '10', zero_inflation = '.1')
  expect_true(is.feature(feat))
  expect_true(feat$get_inter_locus_var())
  par_expr <- feat$get_table()$parameter
  sim <- sapply(1:1000, function(x) {
    locus <- x %% 100;
    eval(parse(text=par_expr))
  })
  expect_true(abs(mean(sim) - theta * 0.9) < .3)
  expect_equal(sum(sim == 5), 0)
  expect_true(sum(sim == 0) > 80)
})


test_that("adding parameter works", {
  feat <- Feature$new()
  expect_equal(feat$add_parameter(5), "5")
  expect_equal(feat$add_parameter("17"), "17")
  expect_equal(feat$add_parameter(par_expr(tau)), "tau")
  expect_equal(length(feat$get_parameters()), 1)
  expect_equal(feat$add_parameter(par_const(5)), "5")

  expect_error(feat$add_parameter(1:2))
  expect_error(feat$add_parameter(NA))
  expect_error(feat$add_parameter(NULL))

  expect_equal(feat$add_parameter(NA, FALSE), NA)
  expect_equal(feat$add_parameter(NULL, FALSE), NA)
})
