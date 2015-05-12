#' @importFrom R6 R6Class
Par_variation <- R6Class("Par_variation", inherit=Parameter,
  private = list(base_par = par_const(NA)),
  public = list(
    initialize = function(parameter, variance) {
      if (is.numeric(parameter) && length(parameter) == 1) {
        expr <- parameter
      } else if (is.character(parameter) && length(parameter) == 1) {
        expr <- parameter
      } else if (is.par(parameter)) {
        idx <- as.character(length(private$parameter) + 1)
        private$base_par <- parameter
        expr <- parameter$get_expression()
      } else {
        stop("Unexpected type of parameter")
      }
      private$expr <- parse(text = paste0('rgamma(1, ', expr, '^2/', variance,
                                          ', ', expr, '/', variance, ')'))
    },
    get_base_par = function() private$base_par
  )
)

is.par_variation <- function(object) inherits(object, "Par_variation")


#' Let the parameter values vary between different loci
#'
#' The function can be used to let the values of a parameter vary between
#' the different loci. When used, the values for the enclosed parameter
#' will follow a gamma distribution with mean of the parameters original
#' value, and the variance specified as argument \code{variance}. This requires
#' that the original value is always positive. When using this, the simulators
#' are called seperately for each locus, which can dramatically increase the
#' time needed to simulate models with many loci.
#'
#' @param par A parameter whichs value will be set variable.
#' @param variance The variance of the gamma distribution, which the values used
#'   for simulation will follow.
#' @export
par_variation <- function(par, variance) {
  Par_variation$new(par, variance)
}

