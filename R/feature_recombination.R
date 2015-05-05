Feature_recombination <- R6Class("Feature_recombination", inherit = Feature,
  private = list(rate = NA),
  public = list(
    initialize = function(rate) {
      private$rate = private$add_parameter(rate)
    },
    get_rate = function() private$rate,
    print = function() {
      cat("Recombination with rate", private$rate, "\n")
    }
  )
)

#' Feature: Recombination
#'
#' Adds intra-locus recombination to a model.
#'
#' The corresponding rate parameter is 4*N0*r, where r is the
#' probability that a recombination event within the locus will
#' occur in one generation. Even when using an infinite sites
#' mutation model, this assumes an finite locus length.
#'
#' @param rate A \code{\link{parameter}} defining the recombination rate
#'   (see above).
#' @return The demographic model with recombination
#' @export
#'
#' @examples
#' # A model with a fixed recombination rate of 5
#' model <- coal_model(c(25,25), 100, 1000) +
#'   feat_recombination(par_const(5))
feat_recombination <- function(rate) {
  Feature_recombination$new(rate)
}

conv_to_ms_arg.Feature_recombination <- function(feature, model) {
  paste0("-r\", format(", feature$get_rate(), ", scientific=FALSE),
         format(locus_length, scientific=FALSE), \"")
}

conv_to_msms_arg.Feature_recombination <- ignore_par
conv_to_seqgen_arg.Feature_recombination <- ignore_par
