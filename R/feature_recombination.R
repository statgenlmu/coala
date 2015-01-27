#' Creates a Recombination feature
#'
#' This function add the assumption to the model that recombination
#' events may occur within each locus. The corresponding parameter
#' equals 4*N0*r, where r is the 
#' probability that a recombination event within the locus will
#' occur in one generation. Even when using an infinite sites
#' mutation model, this assumes an finite locus length which is given
#' by the 'seq.length' parameter of the demographic model.
#' 
#' @inheritParams feat_mutation
#' @param rate A \code{\link{Parameter}} defining the recombination rate 
#'   (see above).
#' @return The demographic model with recombination
#' @export
#'
#' @examples
#' # A model with a fixed recombination rate of 5
#' dm <- dm.createDemographicModel(c(25,25), 100, 1000) + 
#'   feat_recombination(par_const(5))
#' 
#' # A model where the recombination is given by a parameter 'rho'
#' dm <- dm.createDemographicModel(c(25,25), 100) + 
#'   feat_recombination(par_range('rho', 1, 5))
feat_recombination <- function(rate, group = 0, variance = 0) {
  Feature$new("recombination", rate, group = group, variance = variance)
}