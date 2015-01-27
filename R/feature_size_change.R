#' Adds an instantaneous change of the population size of one
#' population to a model.
#'
#' This function changes the effective population size of one
#' population. The change is performed at a given time point
#' ('at.time') and applies to the time interval farther into
#' the past from this point. The population size is set to a
#' fraction of N0, the present day size of population one.
#'
#' If you want to add a slow, continuous change over some time,
#' then use the \link{dm.addGrowth} function.
#'
#' @param dm  The demographic model to which the size change should be added.
#' @param new_size A \code{\link{Parameter}} giving the new size of the population
#'   as fraction of N0.
#' @param population The number of the population whichs size is changed.
#' @param at.time The time point at which the size changes.
#' @return A feature which can be added to a model.
#' @export
#' @examples
#' # A model with one smaller population
#' dm <- dm.createDemographicModel(c(20,37), 88) +
#'   feat_size_change(par_const(.1), 2, at.time='1')
feat_size_change <- function(new_size, population, at.time="0") {
  Feature$new("size.change", new_size, pop_source=population,
              time_point=at.time)
}
