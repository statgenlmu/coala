# Base class for all parameters.
# Contains an expression that can be assigned to some part of a feature.
#' @importFrom R6 R6Class
Parameter <- R6Class('Parameter',  inherit = Base_Object,
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


# Base class for Model Parameters.
# Model parameters have a name, and a value is assigned to a variable of that
# name for each simulation.
#' @importFrom R6 R6Class
Par_Model <- R6Class('Par_Model', inherit=Parameter,
  private = list(name = NA),
  public = list(
    initialize = function(name) {
      if (!(is.character(name) & length(name) == 1))
        stop('The parameter name must be a character')

      super$initialize(parse(text=name))
      private$name <- name
    },
    get_name = function() private$name
  )
)

is.par <- function(par) {
  'Parameter' %in% class(par)
}

is.par_model <- function(par) {
  'Par_Model' %in% class(par)
}

#' Define Model Parameters
#'
#' bla bla bla
#'
#' @param expr An R command. The command will not be evaluted until a simulation
#'  is performed. I can contain other named parameters, but not parameters
#'  created with \code{par_expr}. Make sure that the expression always evaluates
#'  to a suitable parameter value (a single numeric in almost all cases).
#' @describeIn par_expr A parameter whichs value is determined by evalutating an
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
#' @param constant The constant value of the parameter
par_const <- function(constant) {
  Parameter$new(as.expression(constant))
}


Par_Range <- R6Class('Par_Range', inherit = Par_Model,
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
    get_range = function() private$range
  )
)

#' @describeIn par_expr Creates an parameter that can take a range of possible
#'  values. Used for creating model parameters for \code{Jaatha}.
#' @export
#' @param name A string. The name of the parameter.
#' @param lower A numeric. The lower boundary of the parameter's range.
#' @param upper A numeric. The upper boundary of the parameter's range.
par_range <- function(name, lower, upper) {
  Par_Range$new(lower, upper, name)
}

