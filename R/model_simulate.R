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
#' model <- coal_model(c(5,10), 20) +
#'   feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
#'   feat_mutation(par_range('theta', 1, 10)) +
#'   sumstat_jsfs()
#'
#' simulate(model, pars=c(1, 5))
simulate.Coalmodel <- function(object, nsim = 1, seed, pars=NULL, ...) {
  check_par_range(object, pars)

  simprog <- determine_simprog(object)

  if (simprog != "groups") {
    return(get_simprog(simprog)$sim_func(object, pars))
  }

  sum_stats <- list(pars=pars)
  for (group in get_groups(object)) {
    grp_model <- get_group_model(object, group)
    sum_stats_grp <- simulate.Coalmodel(grp_model, pars = pars)

    for (i in seq(along = sum_stats_grp)) {
      if (names(sum_stats_grp)[i] == 'pars') next
      name <- paste(names(sum_stats_grp)[i], group, sep='.')
      sum_stats[[name]] <- sum_stats_grp[[i]]
    }
  }

  sum_stats
}
