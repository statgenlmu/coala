context('Feature Class')

test_that('Creating features works', {
  feat <- Feature$new('abc', 5)
  expect_true(is.feature(feat))
  expect_equal(feat$get_table(), createFeatureTable('abc', '5', NA, NA, NA, 0))
  expect_false(feat$get_inter_locus_var())
  expect_equal(feat$get_group(), 0)
  expect_equal(feat$get_parameters(), list())

  feat <- Feature$new('a', par_expr(2*x), 2, 1, 't5', 1)
  expect_true(is.feature(feat))
  expect_equal(feat$get_table(), createFeatureTable('a', '2 * x', 2, 1, 't5', 1))
  expect_false(feat$get_inter_locus_var())
  expect_equal(feat$get_group(), 1)
  x <- 5
  expect_equal(feat$get_parameters()[["1"]]$eval(), 10)
  
  # Test with two parameters
  feat <- Feature$new('blub', parameter=par_range('a', 5, 7), 
                      time_point=par_range('b', 1, 1.5))
  expect_equal(length(feat$get_parameters()), 2)
  expect_equal(feat$get_table(), createFeatureTable('blub', 'a', NA, NA, 'b', 0))
  
  # Test variance
  feat <- Feature$new('mutation', par_range('theta', 5, 7), variance='10')
  expect_true(is.feature(feat))
  expect_true(feat$get_inter_locus_var())

  par_expr <- feat$get_table()$parameter
  theta = 5
  expect_true(eval(parse(text=par_expr)) != 5)
  sim <- sapply(1:1000, function(x) eval(parse(text=par_expr)))
  expect_true(abs(mean(sim) - theta) < .3)
  
  # Test zero.inflation
  feat <- Feature$new('mutation', par_range('theta', 5, 7), zero_inflation='.1')
  expect_true(is.feature(feat))
  expect_true(feat$get_inter_locus_var())
  
  par_expr <- feat$get_table()$parameter
  dm <- dm.createDemographicModel(5:6, 100)
  locus <- 1; expect_equal(eval(parse(text=par_expr)), 0)
  locus <- 5; expect_equal(eval(parse(text=par_expr)), 0)
  locus <- 10; expect_equal(eval(parse(text=par_expr)), 0)
  locus <- 11; expect_equal(eval(parse(text=par_expr)), 5)
  locus <- 30; expect_equal(eval(parse(text=par_expr)), 5)
  locus <- 72; expect_equal(eval(parse(text=par_expr)), 5)
  
  # Test zero.inflation & variance
  feat <- Feature$new('mutation', parameter = par_range('theta', 5, 7), 
                      variance = '10', zero_inflation = '.1')
  expect_true(is.feature(feat))
  expect_true(feat$get_inter_locus_var())
  par_expr <- feat$get_table()$parameter
  sim <- sapply(1:1000, function(x) { 
    locus <- x %% 100; 
    eval(parse(text=par_expr))
  })
  expect_true(abs(mean(sim) - theta*0.9) < .3)
  expect_equal(sum(sim == 5), 0)
  expect_true(sum(sim == 0) > 80)
})
