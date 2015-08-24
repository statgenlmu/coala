#' @importFrom R6 R6Class
zero_inflation_par_class <- R6Class("zero_inflation_par",
  inherit = variation_par_class,
  private = list(func = "zero_inflate")
)

#' @importFrom stats rbinom
zero_inflate <- function(x, prob) ifelse(rbinom(1, 1, prob), 0, x)


#' @importFrom R6 R6Class
zero_frac_par_class <- R6Class("zero_inflation_par",
  inherit = variation_par_class,
  private = list(func = "zero_frac")
)

create_zero_frac_func <- function(locus_id, locus_number) {
  function(x, frac) ifelse(locus_id <= frac * locus_number, 0, x)
}


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
#'   for each locus if \code{random} is \code{TRUE}. Otherwise, it's the
#'   fixed fraction of loci which will have a parameter value of `0`.
#' @param random Whether the number of loci which are simulated with a value of
#'   `0` should be random (\code{TRUE}) or a fixed fraction (\code{FALSE}). See
#'   \code{prob} for details.
#' @export
par_zero_inflation <- function(par, prob, random = TRUE) {
  if (random) zero_inflation_par_class$new(par, prob)
  else zero_frac_par_class$new(par, prob)
}
