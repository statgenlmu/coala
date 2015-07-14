#' @importFrom R6 R6Class
zero_inflation_par_class <- R6Class("zero_inflation_par",
  inherit = variation_par_class,
  private = list(func = "zero_inflate")
)

#' @importFrom stats rbinom
zero_inflate <- function(x, prob) ifelse(rbinom(1, 1, prob), 0, x)

#' Zero inflation for Parameters
#'
#' This adds a zero inflation to the distribution of a parameter for the
#' different loci. When using this, each locus will be simulated with a
#' parameter value of 0 with probability \code{prob}, or with parameter's
#' original value in the remaining cases.
#  When using this, the simulators
#' are called separately for each locus, which can dramatically increase the
#' time needed to simulate models with many loci.
#'
#' @param par A parameter whose value will be set variable.
#' @param prob The probability that the parameters value will be set to `0`
#'   for each locus.
#' @export
par_zero_inflation <- function(par, prob) {
  zero_inflation_par_class$new(par, prob)
}
