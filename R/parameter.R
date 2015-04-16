# Base class for all parameters.
# Contains an expression that can be assigned to some part of a feature.
#' @importFrom R6 R6Class
Parameter <- R6Class('Parameter',
  private = list(
    expr = NA
  ),
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
  'Parameter' %in% class(par)
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
#' bla bla bla
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
#' @aliases ModelParameters
#' @author Paul Staab
#' @examples
#' par_const(5)         # The parameters value is always 5.
#' par_expr(runif(1))   # Creates an parameter which takes a uniformly
#'                      # distributed value in each simulation.
#' par_range('x', 1, 5) # Creates an parameter with name x with possible values
#'                      # between 1 and 5.
#' par_expr(2*x)        # The parameters value is always equal two twice the
#'                      # value of a different model parameter named 'x'.
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


#' @describeIn par_expr Creates an parameter which's value is specified via the
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


create_par_env <- function(model, parameters, ..., for_cmd = FALSE) {
  par_env <- new.env()

  if (!for_cmd) {
    if (is.null(names(parameters))) {
      par_names <- get_par_names(model)
      if (length(parameters) != length(par_names)) {
        stop("Unexpected number of parameters")
      }
      names(parameters) <- par_names
    }

    for (par in get_parameter(model)) {
      par_env[[par$get_name()]] <- par$generate_value(parameters)
    }
  } else {
    for (par in get_parameter(model)) {
      par_env[[par$get_name()]] <- par$get_name()
    }
  }

  additional_pars <- list(...)
  for (i in seq(along = additional_pars)) {
    par_env[[names(additional_pars)[i]]] <- additional_pars[[i]]
  }

  par_env[["locus_number"]] <- get_locus_number(model)

  par_env
}


escape_par_expr <- function(cmd) {
  if (is.null(names(cmd))) return(cmd)
  for (i in seq(along = cmd)) {
    if (names(cmd)[i] == '') next
    cmd[i] <- paste0('\"', cmd[i], '\"')
  }
  cmd
}

