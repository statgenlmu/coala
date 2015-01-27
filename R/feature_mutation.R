#' Creates a Mutation Feature
#' 
#' This functions adds the assumption to the model that neutral mutations
#' occur in the genomes at a constant rate. The rate is quantified through
#' a parameter usually named theta in population genetics. It equals 4*N0*mu,
#' where N0 is the effective diploid population size of population one at the
#' time of sampling and mu is the neutral mutation rate for an entire locus.
#'
#' @param rate A \code{\link{Parameter}} defining the mutation rate (see above).
#' @param variance Set to a value different from 0 to introduce variation in the
#'                 the parameter value for different loci. The 
#'                 variation follows a gamma distribution with mean equal to
#'                 the value provided as \code{parameter}, and variance as given
#'                 here. Can also be set to a previously 
#'                 created parameter, or an expression based on parameters. 
#' @param group    Group of loci for with this feature is added. 0 means that
#'                 the feature applies to all groups, and 1 is the default group.
#'                 Set to 1 or an greater integer to set this feature only for 
#'                 the corresponding group of loci.
#' @return The feature, which can be added to a model using `+`.
#' @export
#'
#' @examples
#' # A model with a constant scaled mutation rate of 5:
#' dm <- dm.createDemographicModel(c(15,20), 100) + feat_mutation(par_const(5))
#' 
#' # A model with a mutation rate that can be estimated with Jaatha:
#' dm <- dm.createDemographicModel(c(15,20), 100) + 
#'   feat_mutation(par_range('theta', 1, 20))
#'
#' # A model with variable gamma distributed mutation rate
#' dm <- dm.createDemographicModel(c(15,20), 100) + 
#'   feat_mutation(par_range('theta', 1, 20), variance=100)
feat_mutation <- function(rate, group = 0, variance = 0) {
  Feature$new('mutation', parameter=rate,  group=group, variance=variance)
}
