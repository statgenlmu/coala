context('Parameter Class')

test_that('Getting and Setting Expressions works', {
  expect_error(Parameter$new(2 * x))
  expect_error(Parameter$new(2))
  expect_error(Parameter$new('2'))
  basic_par <- Parameter$new(expression(2 * x))
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
  expect_false(is.named_par(basic_par))
})


test_that('par_expr works', {
  basic_par <- par_expr(2 * x)
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
  basic_par <- par_const(2 * x)
  x <- 6
  expect_equal(basic_par$eval(), 10)
})


test_that("par_named works", {
  par <- par_named('theta')
  expect_true(is.par(par))
  expect_true(is.named_par(par))

  expect_equal(par$get_name(), "theta")
  theta <- 5
  expect_equal(par$eval(), 5)

  expect_equal(5, par$generate_value(c(theta = 5)))
  expect_equal(5, par$generate_value(c(x = 2, theta = 5)))
  expect_error(par$generate_value(5))

  expect_true(par$check_value(1))
  expect_true(par$check_value(2))
})


test_that('par_range works', {
  par <- par_range('theta', 1, 2)
  expect_true(is.par(par))
  expect_true(is.named_par(par))

  expect_equal(par$get_name(), "theta")
  theta <- 5
  expect_equal(par$eval(), 5)

  expect_equal(par$get_range(), 1:2)

  expect_true(par$check_value(1))
  expect_true(par$check_value(1.4))
  expect_true(par$check_value(2))
  expect_true(par$check_value(1 - 1e-10))
  expect_true(par$check_value(2 + 1e-10))

  expect_error(par$check_value("1"))
  expect_error(par$check_value(0))
  expect_error(par$check_value(3))
  expect_error(par$check_value(1:2))

  expect_error(par_range('theta', 1:2))
  expect_error(par_range('theta', 2, 1))
})


test_that("Adding an expression par to a model give no error", {
  coal_model(5:6, 10, 100) + par_expr(2 * theta)
  coal_model(5:6, 10, 100) + par_expr(2 * theta) + par_expr(5)
})


test_that('Adding ranged parameters to a model works', {
  model <- coal_model(5:6, 10, 100) + par_range('theta', 1, 2)
  expect_equal(get_parameter_table(model),
               data.frame(name = 'theta', lower.range = 1, upper.range = 2,
                          stringsAsFactors = FALSE))

  model <- coal_model(5:6, 10, 100) +
    par_range('theta', 1, 2) +
    par_range('tau', 5, 6)
  expect_equal(get_parameter_table(model),
               data.frame(name = c('theta','tau'),
                          lower.range = c(1, 5),
                          upper.range = c(2, 6),
                          stringsAsFactors = FALSE))
})


test_that("Creation of parameter enviroment works", {
  # With named parameters
  model <- coal_model() + par_named("x")
  par_envir <- create_par_env(model, c(x = 5))
  expect_equal(par_envir[['x']], 5)
  par_envir <- create_par_env(model, 5)
  expect_equal(par_envir[['x']], 5)
  par_envir <- create_par_env(model, c(y = 2, x = 5))
  expect_equal(par_envir[['x']], 5)
  expect_error(create_par_env(model, numeric()))
  expect_error(create_par_env(model, 1:2))
  expect_error(create_par_env(model, c(y = 2)))

  # Without parameters
  par_envir <- create_par_env(coal_model(), numeric(0))

  # With ranged parameters (not really needed)
  par_envir <- create_par_env(model_theta_tau(), c(1, 5))
  expect_equal(par_envir[['tau']], 1)
  expect_equal(par_envir[['theta']], 5)

  par_envir <- create_par_env(model_theta_tau(), c(theta = 5, tau = 1))
  expect_equal(par_envir[['tau']], 1)
  expect_equal(par_envir[['theta']], 5)

  # Additional options
  par_envir <- create_par_env(model_theta_tau(), c(1, 5), locus = 17)
  expect_equal(par_envir[['locus']], 17)

  par_envir <- create_par_env(model_theta_tau(), c(1, 5),
                              locus = 23, seed = 115)
  expect_equal(par_envir[['locus']], 23)
  expect_equal(par_envir[['seed']], 115)


  # For cmd printing
  par_envir <- create_par_env(model_theta_tau(), for_cmd = TRUE)
  expect_equal(par_envir[["theta"]], "theta")
  expect_equal(par_envir[["tau"]], "tau")
})


