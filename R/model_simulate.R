#' Simulates data according to a demographic model
#'
#' @param dm The demographic model according to which the simulations should be done
#' @param parameters A vector of parameters which should be used for the simulations.
#'           If a matrix is given, a simulation for each row of the matrix will be performed
#' @return A matrix where each row is the vector of summary statistics for
#'         the parameters in the same row of the "parameter" matrix
#' @export
#'
#' @examples
#' model <- CoalModel(c(5,10), 20) +
#'   feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
#'   feat_mutation(par_range('theta', 1, 10)) +
#'   sumstat_jsfs()
#'
#' simulate(model, c(1, 5))
simulate.CoalModel <- function(dm, parameters) {
  stopifnot(is.model(dm))
  checkParInRange(dm, parameters)

  if (!dm$finalized) dm = dm.finalize(dm)

  if (dm$currentSimProg != "groups") {
    return(getSimProgram(dm$currentSimProg)$sim_func(dm, parameters))
  }

  sum_stats <- list(pars=parameters)
  for (group in get_groups(dm)) {
    dm.grp <- dm$options$grp.models[[as.character(group)]]
    sum_stats.grp <- getSimProgram(dm.grp$currentSimProg)$sim_func(dm.grp, parameters)
    for (i in seq(along = sum_stats.grp)) {
      if (names(sum_stats.grp)[i] == 'pars') next()
      name <- paste(names(sum_stats.grp)[i], group, sep='.')
      sum_stats[[name]] <- sum_stats.grp[[i]]
    }
  }

  sum_stats
}
