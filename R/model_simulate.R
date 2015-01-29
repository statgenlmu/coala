#' Simulates data according to a demographic model
#'
#' @param dm The demographic model according to which the simulations should be done
#' @param parameters A vector of parameters which should be used for the simulations.
#'           If a matrix is given, a simulation for each row of the matrix will be performed
#' @return A matrix where each row is the vector of summary statistics for
#'         the parameters in the same row of the "parameter" matrix
#' @export
#' @importFrom stats simulate
#'
#' @examples
#' model <- CoalModel(c(5,10), 20) +
#'   feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
#'   feat_mutation(par_range('theta', 1, 10)) +
#'   sumstat_jsfs()
#'
#' simulate(model, pars=c(1, 5))
simulate.CoalModel <- function(object, nsim = 1, seed, pars=NULL, ...) {
  stopifnot(!is.null(pars))
  checkParInRange(object, pars)

  if (!object$finalized) object = dm.finalize(object)

  if (object$currentSimProg != "groups") {
    return(getSimProgram(object$currentSimProg)$sim_func(object, pars))
  }

  sum_stats <- list(pars=pars)
  for (group in get_groups(object)) {
    object.grp <- object$options$grp.models[[as.character(group)]]
    sum_stats.grp <- getSimProgram(object.grp$currentSimProg)$sim_func(object.grp, pars)
    for (i in seq(along = sum_stats.grp)) {
      if (names(sum_stats.grp)[i] == 'pars') next()
      name <- paste(names(sum_stats.grp)[i], group, sep='.')
      sum_stats[[name]] <- sum_stats.grp[[i]]
    }
  }

  sum_stats
}
