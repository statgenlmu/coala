recombination_class <- R6Class("recombination", inherit = feature_class,
  public = list(
    initialize = function(rate) {
      private$rate <- private$add_parameter(rate)
    },
    print = function() {
      cat("Recombination with rate", print_par(private$rate), "\n")
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
#'   feat_recombination(5)
feat_recombination <- function(rate) {
  recombination_class$new(rate)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.recombination <- function(feature, model) {
  paste0("-r', ", feature$get_rate(), ", par(locus_length), '")
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.recombination <- conv_to_ms_arg.recombination

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.recombination <- conv_to_ms_arg.recombination

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.recombination <- ignore_par
