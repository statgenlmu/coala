#' Creates a feature that represents the sampling from one population
#'
#' @param size The number of individuals that are sampled.
#' @param population The population from with the indidivuals are sampled
#' @param time_point The time at which the sample is taken.
#' @return The feature, which can be added to a model using `+`.
#' @export
feat_unphased <- function(ploidy, samples_per_ind=ploidy) {
  Feature$new('unphased', par_const(NA)) +
    Feature$new('ploidy', par_const(ploidy)) +
    Feature$new('samples_per_ind', par_const(samples_per_ind))
}
