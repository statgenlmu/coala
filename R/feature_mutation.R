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
#' @param model    The mutation model you want to use. Can be IFS (default),
#'                 HKY, F84 or GTR.
#' @param tstv.ratio The ratio of transitions to transversions. The default is
#'                   0.5, which means that all amino acid substitutions are
#'                   equally likely. In this case, the HKY model is identical to
#'                   the Felsenstein 81 model.
#' @param base.frequencies The equilibrium frequencies of the four bases.
#'                   Must be a numeric vector of length four.
#'                   Order is A, C, G, T.
#' @param gtr.rates  The rates for the amino acid substitutions. Must be a
#'                   numeric vector of length six. Order: A->C, A->G, A->T, C->G, C->T, G->T.
#' @return The feature, which can be added to a model using `+`.
#' @export
#'
#' @examples
#' # A model with a constant scaled mutation rate of 5:
#' dm <- coal_model(c(15,20), 100) + feat_mutation(par_const(5))
#'
#' # A model with a mutation rate that can be estimated with Jaatha:
#' dm <- coal_model(c(15,20), 100) +
#'   feat_mutation(par_range('theta', 1, 20))
#'
#' # A model with variable gamma distributed mutation rate
#' dm <- coal_model(c(15,20), 100) +
#'   feat_mutation(par_range('theta', 1, 20), variance=100)
feat_mutation <- function(rate, variance = 0, model='IFS',
                          base_frequencies, tstv_ratio, gtr_rates) {

  container <- coal_model() + Feature$new('mutation', parameter=rate,
                                          variance=variance)

  # Add the mutation model
  if (model != 'IFS') {
    if (!model %in% sg_mutation_models)
      stop("Possible mutation models: ",
           paste(sg_mutation_models, collapse=" "))
    container <- container + Feature$new("mutation_model", par_const(model))
  }

  if ( !missing(tstv_ratio) ) {
    if (!model %in% c("HKY", "F84"))
      stop("This mutation model does not support a ts/tv ratio")
    container <- container + Feature$new("tstv_ratio", tstv_ratio)
  }

  if ( !missing(base_frequencies) ) {
    stopifnot(length(base_frequencies) == 4)
    stopifnot(sum(base_frequencies) == 1)
    if (!model %in% c("HKY", "F84"))
      stop("This mutation model does not support base frequencies")

    container <- container + Feature$new("base_freq_A", base_frequencies[1])
    container <- container + Feature$new("base_freq_C", base_frequencies[2])
    container <- container + Feature$new("base_freq_G", base_frequencies[3])
    container <- container + Feature$new("base_freq_T", base_frequencies[4])
  }

  if (!missing(gtr_rates)) {
    stopifnot(length(gtr_rates) == 6)
    if (!model %in% c("GTR"))
      stop("You can specify gtr_rates only with the GTR model")

    container <- container + Feature$new("gtr_rate_1", gtr_rates[1])
    container <- container + Feature$new("gtr_rate_2", gtr_rates[2])
    container <- container + Feature$new("gtr_rate_3", gtr_rates[3])
    container <- container + Feature$new("gtr_rate_4", gtr_rates[4])
    container <- container + Feature$new("gtr_rate_5", gtr_rates[5])
    container <- container + Feature$new("gtr_rate_6", gtr_rates[6])
  }

  container
}

#-------------------------------------------------------------------
# dm.setMutationModel
#-------------------------------------------------------------------
# Defines what mutation model is used for simulations
#
# As default, we simulate mutation using the Infinite Sites Model.
# Using the function, you can change it either to the Hasegawa, Kishino and
# Yano (HKY), to the Felsenstein and Churchill 96 (F84) or to the Generalised
# time reversible (GTR) model.
# This requires that seq-gen is installed on our system.
#
# The HKY and F84 models use the the arguments 'base.frequencies' and
# 'tstv.ratio'. The GTR model uses 'gtr.rates'.


# dm.addMutationRateHeterogenity <-
#   function(dm, min.alpha, max.alpha, parameter="alpha", categories.number) {
#
#     dm <- addFeature(dm, "gamma.rate", parameter, min.alpha,
#                      max.alpha, NA, NA, NA)
#
#     if (!missing(categories.number)) {
#       dm <- addFeature(dm, "gamma.categories", parameter=categories.number)
#     }
#
#     return(dm)
#   }


# Set the mutation rates for trios
# @param middle_rate The mutation rate used for the middle locus
# @param outer_rate The mutation rate for the two outer loci
# @export
# dm.setTrioMutationRates <- function(dm, middle_rate, outer_rate, group = 0) {
#   dm <- addFeature(dm, 'mutation', parameter = middle_rate, group = group)
#   dm <- addFeature(dm, 'mutation_outer', parameter = outer_rate,
#                    group = group)
# }
