#' Creates a feature that represents the sampling from one population
#' 
#' @param size The number of individuals that are sampled.
#' @param population The population from with the indidivuals are sampled
#' @param time_point The time at which the sample is taken.
#' @return The feature, which can be added to a model using `+`.
#' @export
feat_sample <- function(size, population, time_point='0') {
  if (time_point != '0')
    stop("Samples at time different from 0 at not supported at the moment")
  Feature$new('sample', size, pop_source=population, time_point=time_point)
}