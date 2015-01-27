context('Parameter Class')

test_that('Getting and Setting Expressions works', {
  expect_error(Parameter$new(2*x))
  expect_error(Parameter$new(2))
  expect_error(Parameter$new('2'))
  basic_par <- Parameter$new(expression(2*x))
  x <- 5
  expect_equal(basic_par$eval(), 10)
  
  test_env <- new.env()
  test_env[['x']] <- 6
  expect_equal(basic_par$eval(envir = test_env), 12)
  expect_equal(basic_par$eval(), 10)
  
  expr <- basic_par$get_expression()
  expect_is(expr, 'expression')
  expect_equal(eval(expr), 10)
  
  expect_true(is.par(basic_par))
  expect_false(is.par_model(basic_par))
})


test_that('par_expr works', {
  basic_par <- par_expr(2*x)
  expect_is(basic_par, 'Parameter')
  x <- 5
  expect_equal(basic_par$eval(), 10)
  x <- 6
  expect_equal(basic_par$eval(), 12)
  
  basic_par <- par_expr(sqrt(y))
  y <- 4
  expect_equal(basic_par$eval(), 2)
})


test_that('par_const works', {
  x <- 5
  basic_par <- par_const(2*x)
  x <- 6
  expect_equal(basic_par$eval(), 10)
})



test_that('par_range works', {
  par <- par_range('theta', 1, 2)
  expect_true(is.par(par))
  expect_true(is.par_model(par))
  
  expect_equal(par$get_name(), "theta")
  expect_equal(par$get_range(), 1:2)
  theta <- 5
  expect_equal(par$eval(), 5)
  
  expect_error(par_range('theta', 1:2))
  expect_error(par_range('theta', 2, 1))
})


test_that("Adding an expression par to a model give no error", {
  dm.createDemographicModel(5:6, 10, 100) + par_expr(2*theta)
  dm.createDemographicModel(5:6, 10, 100) + par_expr(2*theta) + par_expr(5)
})


test_that('Adding ranged parameters to a model works', {
  model <- dm.createDemographicModel(5:6, 10, 100) + par_range('theta', 1, 2)
  expect_equal(model@parameters, 
               data.frame(name='theta', lower.range=1, upper.range=2, 
                          stringsAsFactors = FALSE))
  
  model <- dm.createDemographicModel(5:6, 10, 100) + 
    par_range('theta', 1, 2) +
    par_range('tau', 5, 6)
  expect_equal(model@parameters, 
               data.frame(name=c('theta','tau'), 
                          lower.range=c(1, 5), 
                          upper.range=c(2, 6), 
                          stringsAsFactors = FALSE))
})