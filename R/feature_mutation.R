#' @importFrom R6 R6Class
mutation_class <- R6Class("mutation", inherit = feature_class,
  private = list(
    model = NA,
    fixed = FALSE,
    base_frequencies = NA,
    tstv_ratio = NA,
    gtr_rates = NA
  ),
  public = list(
    initialize = function(rate, model, base_frequencies,
                          tstv_ratio, gtr_rates, fixed) {
      private$rate <- private$add_parameter(rate, add_par = FALSE)

      assert_that(length(model) == 1)
      assert_that(any(model == c("IFS", "HKY", "GTR")))
      private$model <- model

      assert_that(length(fixed) == 1)
      assert_that(is.logical(fixed))
      private$fixed <- fixed

      if (model == "HKY") {
        if (is.na(tstv_ratio)) {
          stop("You need to specify tstv_ratio for the HKY mutation model")
        }
        assert_that(all(!is.na(base_frequencies)))
        assert_that(is.numeric(tstv_ratio))
        assert_that(length(tstv_ratio) == 1)
        private$tstv_ratio <- tstv_ratio

        if (any(is.na(base_frequencies))) {
          stop("missing base_frequencies for the HKY mutation model")
        }
        assert_that(all(!is.na(base_frequencies)))
        assert_that(is.numeric(base_frequencies))
        assert_that(length(base_frequencies) == 4)
        assert_that(sum(base_frequencies) == 1)
        private$base_frequencies <- base_frequencies
      }

      else if (model == "GTR") {
        if (any(is.na(gtr_rates))) {
          stop("You need to specify gtr_rates for the GTR mutation model")
        }
        assert_that(is.numeric(gtr_rates))
        assert_that(length(gtr_rates) == 6)
        private$gtr_rates <- gtr_rates
      }
    },
    get_model = function() private$model,
    get_base_frequencies = function() private$base_frequencies,
    get_tstv_ratio = function() private$tstv_ratio,
    get_gtr_rates = function() private$gtr_rates,
    get_fixed = function() private$fixed,
    print = function() {
      cat("Mutations with rate", print_par(private$rate),
          "following a", private$model, "mutation model\n")
    }
  )
)

#' Feature: Mutation
#'
#' This functions adds the assumption to the model that neutral mutations
#' occur in the genomes at a constant rate. The rate is quantified through
#' a parameter usually named theta in population genetics. It equals 4*N0*mu,
#' where N0 is the effective diploid population size of population one at the
#' time of sampling and mu is the neutral mutation rate for an entire locus.
#'
#' @param rate A \code{\link{parameter}} defining the mutation rate.
#' @param fixed_number If set to \code{TRUE}, the number of mutations on each
#'   locus will always be exactly equal to the rate, rather than happening with
#'   a rate along the ancestral tree.
#' @param model The mutation model you want to use.
#'   Can be either 'IFS' (default), 'HKY' or 'GTR'. Refer to the mutation model
#'   section for detailed information.
#' @param tstv_ratio The ratio of transitions to transversions used in the 'HKY'
#'   muation model.
#' @param base_frequencies The equilibrium frequencies of the four bases used in
#'   the 'HKY' mutation model. Must be a numeric vector of length four, with the
#'   values for A, C, G and T, in that order.
#' @param gtr_rates The rates for the six amino acid substitutions used in the
#'   'GTR' model. Must be a numeric vector of length six.
#'   Order: A<->C, A<->G, A<->T, C<->G, C<->T, G<->T.
#' @return The feature, which can be added to a model using `+`.
#' @export
#'
#' @section Mutation Models:
#' The Hasegawa, Kishino and Yano (HKY) model (Hasegawa et al., 1985) allows
#' for a different rate of transitions and transversions (tstv_ratio)
#' and unequal
#' frequencies of the four nucleotides (base_frequencies).
#'
#' The general reversible process (GTR) model (e.g. Yang, 1994) is more general
#' than the HKY model and allows to define the rates for each
#' type of substitution. The rates are assumed to be symmetric
#' (e.g., the rate for T to G is equal to the one for G to T).
#'
#' @examples
#' # A model with a constant mutation rate of 5:
#' model <- coal_model(10, 1) + feat_mutation(rate = 5)
#'
#' # A model with a mutation rate that can be estimated with Jaatha:
#' model <- coal_model(c(15,20), 100) +
#'   feat_mutation(par_range('theta', 1, 20))
feat_mutation <- function(rate,
                          model = "IFS",
                          base_frequencies = NA,
                          tstv_ratio = NA,
                          gtr_rates = NA,
                          fixed_number = FALSE) {

  mutation_class$new(rate, model, base_frequencies,
                     tstv_ratio, gtr_rates, fixed_number)
}

#-------------------------------------------------------------------
# model.setMutationModel
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


# model.addMutationRateHeterogenity <-
#   function(model, min.alpha, max.alpha, parameter="alpha", categories.number) {
#
#     model <- addFeature(model, "gamma.rate", parameter, min.alpha,
#                      max.alpha, NA, NA, NA)
#
#     if (!missing(categories.number)) {
#       model <- addFeature(model, "gamma.categories", parameter=categories.number)
#     }
#
#     return(model)
#   }


# Set the mutation rates for trios
# @param middle_rate The mutation rate used for the middle locus
# @param outer_rate The mutation rate for the two outer loci
# @export
# model.setTrioMutationRates <- function(model, middle_rate, outer_rate, group = 0) {
#   model <- addFeature(model, 'mutation', parameter = middle_rate, group = group)
#   model <- addFeature(model, 'mutation_outer', parameter = outer_rate,
#                    group = group)
# }

is_feat_mutation <- function(feat) any("mutation" == class(feat))

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.mutation <- function(feature, model) {
  if (feature$get_model() != "IFS") stop("Unsupported mutation model")
  if (feature$get_fixed()) {
    paste0("-s', par(", feature$get_rate(), "), '")
  } else {
    paste0("-t', par(", feature$get_rate(), "), '")
  }
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.mutation <- conv_to_ms_arg.mutation

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.mutation <- function(feature, model) {
  if (feature$get_fixed()) {
    stop("scrm does not support simulating a fixed number of mutations")
  }
  conv_to_ms_arg.mutation(feature, model)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.mutation <- function(feature, model) {
  if (feature$get_model() == "IFS") {
    stop("seq-gen can not simulate an IFS model")
  }
  if (feature$get_fixed()) {
    stop("seq-gen can not simulate a fixed number of mutations")
  }
  if (feature$get_model() == "GTR") {
    rates <- paste("-r", paste(feature$get_gtr_rates(), collapse = " "))
  } else {
    rates <- paste("-f", paste(feature$get_base_frequencies(), collapse = " "),
                   "-t", feature$get_tstv_ratio())
  }
  paste0("-m", feature$get_model(), " ",
         rates, " ",
         "-l', locus_length, '",
         "-s', par(", feature$get_rate(), " / locus_length), '",
         "-p', locus_length + 1, '",
         "-z', par(seed), '-q")
}
