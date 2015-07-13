# Base class for all parameters.
# Contains an expression that can be assigned to some part of a feature.
#' @importFrom R6 R6Class
Parameter <- R6Class('Parameter',
  private = list(expr = NA),
  public = list(
    initialize = function(expr) {
      if (!is.expression(expr)) stop("No expression provided: ", expr,
                                     " is of type ", is(expr))
      private$expr <- expr
    },
    eval = function(envir = parent.frame()) {
      eval(private$expr, envir = envir)
    },
    get_expression = function() private$expr
  )
)


is.par <- function(par) {
  any('Parameter' == class(par))
}



# Base class for Model Parameters.
# Model parameters have a name, and a value is assigned to a variable of that
# name for each simulation.
#' @importFrom R6 R6Class
Par_Named <- R6Class('Par_Named', inherit=Parameter,
  private = list(name = NA),
  public = list(
    initialize = function(name) {
      if (!(is.character(name) & length(name) == 1))
        stop('The parameter name must be a character')

      super$initialize(parse(text=name))
      private$name <- name
    },
    get_name = function() private$name,
    print = function() {
      cat(private$name, ": Named parameter")
    },
    check_value = function(value) TRUE,
    generate_value = function(pars) {
      if (is.null(names(pars)) || !any(private$name == names(pars))) {
        stop("No value for parameter ", private$name, " found")
      }
      value <- pars[[private$name]]
      self$check_value(value)
      value
    }
  )
)


is.named_par <- function(par) {
  'Par_Named' %in% class(par)
}


#' Define Model Parameters
#'
#' This functions allow to add parameters to a model. Parameters can either
#' be used in features, or added directly to a model using the plus operator.
#' The value of parameters can be specified in the simulation command
#' (for \code{par_named} and \code{par_range}), sampled from a prior
#' distribution (\code{par_prior}) or can be derived from other parameters
#' (\code{par_expr}).
#'
#' @param expr An R expression.
#'  This allows to define a parameter using an R expression.
#'  It can contain other named parameters (e.g. \code{2 * a} will create an
#'  parameter that is twice the value of an existing parameter \code{a}).
#'  Make sure that the expression always evaluates
#'  to a valid parameter value (a single numeric in almost all cases).
#' @describeIn par_expr Creates a parameter with value determined by evaluating an
#'  expression.
#' @export
#' @aliases parameter
#' @author Paul Staab
#' @examples
#' par_const(5)
#' par_named("x")
#' par_prior("y", rnorm(1))
#' par_range("z", 1, 5)
#' par_expr(2*x + y * z)
par_expr <- function(expr) {
  Parameter$new(as.expression(substitute(expr)))
}


#' @describeIn par_expr Creates an parameter that is equal to a fixed value.
#'   Different to par_expr, the value is evaluated on parameter creation.
#' @export
#' @param constant An R expression.
#'   The constant value of the parameter.
#'   Different to \code{expr}, the expression is evaluated immediately and
#'   can not depend on other named parameters.
par_const <- function(constant) {
  Parameter$new(as.expression(constant))
}


#' @describeIn par_expr Creates an parameter whose value is specified via the
#'   \code{pars} argument in \code{\link{simulate.Coalmodel}}.
#' @export
#' @param name Character. The name of the parameter. Must be unique in a model.
par_named <- function(name) {
  Par_Named$new(name)
}


Par_Range <- R6Class('Par_Range', inherit = Par_Named,
  private = list(range = NA),
  public = list(
    initialize = function(lower, upper, name) {
      stopifnot(all(is.numeric(c(lower, upper))))
      stopifnot(length(lower) == 1)
      stopifnot(length(upper) == 1)
      stopifnot(lower < upper)

      super$initialize(name)
      private$range <- c(lower, upper)
    },
    get_range = function() private$range,
    print = function() {
      cat(private$name, ": range between", private$range[1],
                           'and', private$range[2], "\n")
    },
    check_value = function(value) {
      if ((!is.numeric(value)) || length(value) != 1) {
        stop("The value for ", private$name, " must be a single numeric value")
      }
      if (value < private$range[1] - 1e-5 || value > private$range[2] + 1e-5) {
        stop("The value for ", private$name, " is not in the given range.")
      }
      TRUE
    }
  )
)


is.ranged_par <- function(par) {
  'Par_Range' %in% class(par)
}


#' @describeIn par_expr Creates an parameter that can take a range of possible
#'  values.
#'  Similar to \code{\link{par_named}}, the value of the parameter
#'  used in a simulation is set via the \code{pars} argument.
#'  This is primarily intended for creating model parameters for
#'  \pkg{jaatha}.
#'
#' @export
#' @param lower A numeric. The lower boundary of the parameter's range.
#' @param upper A numeric. The upper boundary of the parameter's range.
par_range <- function(name, lower, upper) {
  Par_Range$new(lower, upper, name)
}


par_eval_func <- function(x) format(x, scientific = FALSE)
par_print_func <- function(x) substitute(x)

create_par_env <- function(model, parameters, ..., for_cmd = FALSE) {
  par_env <- new.env()

  if (!for_cmd) {
    par_env[['par']] <- par_eval_func

    for (par in get_parameter(model)) {
      par_env[[par$get_name()]] <- par$generate_value(parameters)
    }
  } else {
    par_env[['par']] <- par_print_func
  }

  additional_pars <- list(...)
  for (i in seq(along = additional_pars)) {
    par_env[[names(additional_pars)[i]]] <- additional_pars[[i]]
  }

  par_env
}


prepare_pars <- function(pars, model) {
  assert_that(is.numeric(pars))

  # Try to add names for non-prior parameters if missing
  if (is.null(names(pars))) {
    par_names <- get_par_names(model, without_priors = TRUE)
    if (length(pars) != length(par_names)) {
      stop("Unexpected number of parameters")
    }
    names(pars) <- par_names
  }

  # Sample from priors and return
  c(pars, sample_par_priors(model))
}


print_par <- function(par) paste0("`", substr(par, 5, nchar(par) - 1), "`")
