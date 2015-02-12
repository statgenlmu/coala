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
#' dm <- CoalModel(c(15,20), 100) + feat_mutation(par_const(5))
#'
#' # A model with a mutation rate that can be estimated with Jaatha:
#' dm <- CoalModel(c(15,20), 100) +
#'   feat_mutation(par_range('theta', 1, 20))
#'
#' # A model with variable gamma distributed mutation rate
#' dm <- CoalModel(c(15,20), 100) +
#'   feat_mutation(par_range('theta', 1, 20), variance=100)
feat_mutation <- function(rate, group = 0, variance = 0, model='IFS',
                          base_frequencies, tstv_ratio, gtr_rates) {
  feat <- Feature$new('mutation', parameter=rate,  group=group, variance=variance)

  # Add the mutation model
  if (model != 'IFS') {
    if (! model %in% sg.mutation.models)
      stop("Possible mutation models: ", paste(sg.mutation.models, collapse=" "))
    feat$add_feature(Feature$new("mutation_model", par_const(model)))
  }

  if ( !missing(tstv_ratio) ) {
    if (!model %in% c("HKY", "F84"))
      stop("This mutation model does not support a ts/tv ratio")
    feat$add_feature(Feature$new("tstv_ratio", tstv_ratio))
  }

  if ( !missing(base_frequencies) ) {
    stopifnot(length(base_frequencies) == 4)
    stopifnot(sum(base_frequencies) == 1)
    if (!model %in% c("HKY", "F84"))
      stop("This mutation model does not support base frequencies")

    feat$add_feature(Feature$new("base_freq_A", base_frequencies[1]))
    feat$add_feature(Feature$new("base_freq_C", base_frequencies[2]))
    feat$add_feature(Feature$new("base_freq_G", base_frequencies[3]))
    feat$add_feature(Feature$new("base_freq_T", base_frequencies[4]))
  }

  if (!missing(gtr_rates)) {
    stopifnot(length(gtr_rates) == 6)
    if (!model %in% c("GTR"))
      stop("You can specify gtr_rates only with the GTR model")

    feat$add_feature(Feature$new("gtr_rate_1", gtr_rates[1]))
    feat$add_feature(Feature$new("gtr_rate_2", gtr_rates[2]))
    feat$add_feature(Feature$new("gtr_rate_3", gtr_rates[3]))
    feat$add_feature(Feature$new("gtr_rate_4", gtr_rates[4]))
    feat$add_feature(Feature$new("gtr_rate_5", gtr_rates[5]))
    feat$add_feature(Feature$new("gtr_rate_6", gtr_rates[6]))
  }

  feat
}

#-------------------------------------------------------------------
# dm.setMutationModel
#-------------------------------------------------------------------
# Defines what mutation model is used for simulations
#
# As default, we simulate mutation using the Infinite Sites Model.
# Using the function, you can change it either to the Hasegawa, Kishino and
# Yano (HKY), to the Felsenstein and Churchill 96 (F84) or to the Generalised
# time reversible (GTR) model. This requires that seq-gen is installed on our system.
#
# The HKY and F84 models use the the arguments 'base.frequencies' and
# 'tstv.ratio'. The GTR model uses 'gtr.rates'.



#-------------------------------------------------------------------
# dm.addMutationRateHeterogenity
#-------------------------------------------------------------------
#' Allows the mutation rate on different sites within one locus to
#' vary according to a Gamma Distribution.
#'
#' This function adds a Gamma distributed rate heterogeneity as implemented
#' in 'seq-gen' to the model.
#'
#' "The [...] model of rate heterogeneity assigns different rates to different
#' sites according to a gamma distribution (Yang, 1993). The distribution is scaled
#' such that the mean rate for all the sites is 1 but
#' the user must supply a parameter which describes its shape. A low value for this
#' parameter (<1.0) simulates a large degree of site-specific rate heterogeneity
#' and as this value increases the simulated data becomes more rate-homogeneous.
#' This can be performed as a continuous model, i.e. every site has a different
#' rate sampled from the gamma distribution of the given shape, or as a discrete
#' model, i.e. each site falls into one of N rate categories approximating the
#' gamma distribution. For a review of site-specific rate heterogeneity and its
#' implications for phylogenetic analyses, see Yang (1996)."
#' [From the seq-gen homepage http://bioweb2.pasteur.fr/docs/seq-gen ]
#'
#' The Parameter in this text will be referred to as 'alpha'. Simulation a model
#' with rate heterogeneity requires that 'seq-gen' is installed on your system.
#'
#' @param dm  The demographic model to which the rate heterogeneity should be added.
#' @param min.alpha  If you want to estimate the rate heterogeneity, this will be
#'            used as the smallest possible value.
#' @param max.alpha  Same as min.growth.rate, but the largest possible value.
#' @param categories.number If this is set, a fixed number of categories will be
#'            used to model the gamma distribution instead of drawing every parameter
#'            seperately (see text).
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "alpha" will use
#'            a parameter with name alpha that you have previously
#'            created. You can also use R expression here, i.e. "2*alpha"
#'            or "5*M+2*alpha" (if M is another parameter) will also
#'            work (also the latter does not make much sense).
#' @return    The demographic model with a size change.
dm.addMutationRateHeterogenity <-
  function(dm, min.alpha, max.alpha, parameter="alpha", categories.number) {

    dm <- addFeature(dm, "gamma.rate", parameter, min.alpha, max.alpha, NA, NA, NA)

    if (!missing(categories.number)) {
      dm <- addFeature(dm, "gamma.categories", parameter=categories.number)
    }

    return(dm)
  }


#' Set the mutation rates for trios
#' @param middle_rate The mutation rate used for the middle locus
#' @param outer_rate The mutation rate for the two outer loci
#' @export
dm.setTrioMutationRates <- function(dm, middle_rate, outer_rate, group = 0) {
  dm <- addFeature(dm, 'mutation', parameter = middle_rate, group = group)
  dm <- addFeature(dm, 'mutation_outer', parameter = outer_rate, group = group)
}
