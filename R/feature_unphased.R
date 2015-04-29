#' Creates a feature that represents the sampling from one population
#'
#' @param size The number of individuals that are sampled.
#' @param population The population from with the indidivuals are sampled
#' @param time_point The time at which the sample is taken.
#' @return The feature, which can be added to a model using `+`.
#' @export
feat_unphased <- function(ploidy, samples_per_ind=ploidy) {
  if (samples_per_ind > ploidy) {
    stop("samples_per_ind can not be larger than the ploidy.")
  }

  coal_model() +
    Feature$new('unphased', par_const(NA)) +
    Feature$new('ploidy', par_const(ploidy)) +
    Feature$new('samples_per_ind', par_const(samples_per_ind))
}


get_ploidy <- function(model) {
  feat_ploidy <- search_feature(model, type="ploidy")
  if (nrow(feat_ploidy) == 0) return(1L)
  assert_that(nrow(feat_ploidy) == 1)
  as.integer(feat_ploidy$parameter)
}


get_samples_per_ind <- function(model) {
  feat_spi <- search_feature(model, type="samples_per_ind")
  if (nrow(feat_spi) == 0) return(1L)
  assert_that(nrow(feat_spi) == 1)
  as.integer(feat_spi$parameter)
}


is_unphased <- function(model) {
  feat_spi <- search_feature(model, type="unphased")
  if (nrow(feat_spi) == 0) return(FALSE)
  TRUE
}
