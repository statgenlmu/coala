#' Simulates data according to a demographic model
#'
#' @param object The coalescent model to be simulated
#' @param nsim currently unused
#' @param seed A random seed that is set before simulation.
#' @param ... currently unused
#' @param pars Values for parameters specified with \code{\link{par_named}} or
#'   \code{\link{par_range}}. Should be a named numeric vector.
#' @return A list of summary statistics.
#' @export
#' @importFrom stats simulate
#' @importFrom assertthat assert_that
#'
#' @examples
#' model <- coal_model(c(5,10), 20) +
#'   feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
#'   feat_mutation(par_range('theta', 1, 10)) +
#'   sumstat_jsfs()
#'
#' simulate(model, pars=c(1, 5))
simulate.coalmodel <- function(object, nsim = 1, seed, ..., pars = numeric(0)) {
  if (!missing(seed)) set.seed(seed)
  simprog <- select_simprog(object)
  if (is.null(simprog)) stop("No simulator found")

  results <- list()
  length(results) <- nsim
  for (i in seq(len = nsim)) {
    current_pars <- prepare_pars(pars, object)
    results[[i]] <- simprog$simulate(object, current_pars)
  }

  if (nsim == 1) return(results[[1]])
  results
}
