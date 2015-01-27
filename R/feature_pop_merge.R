#' Adds a speciation event to a demographic model
#'
#' You can use this function the create a new population that splits 
#' of from an existing population at a given time in the past. The time
#' can be given as parameter or as an expression based on previously
#' generated time points.
#'
#' Time in measured in Number of 4N0 generations in the past,
#' where N0 is the size of population 1 at time 0.
#'
#' @param in.pop The population in which the spilt
#'            occurs. See above for more information.
#' @param to.pop The newly created population.
#' @param time.point  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "tau" will use
#'            an parameter with name tau that you have previously 
#'            created. You can also use R expression here, i.e. "2*tau"
#'            or "5*M+2*tau" (if M is another parameter) will also
#'            work (also this does not make much sense).
#' @return    The demographic model with a split.
#' @export
#' @examples
#' dm <- dm.createDemographicModel(c(25,25), 100) + 
#'   feat_pop_merge(par_range('tau', 0.01, 5), 2, 1) +
#'   feat_mutation(par_range('theta', 1, 20))
feat_pop_merge <- function(time, pop_source, pop_target) {
  Feature$new('pop_merge', parameter=par_const(NA), pop_source=pop_source, 
              pop_sink=pop_target, time_point=time)
}