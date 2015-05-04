Feature_outgroup <- R6Class("Feature_outgroup", inherit = Feature,
  public = list(
    initialize = function(population) {
      private$set_population(population)
    },
    print = function() {
      cat("Outgroup: Population", private$population, "\n")
    }
  )
)

#' Adds an outgroup to a demographic model
#'
#' This feature declares an existing population as outgroup. Outgroups are used
#' to determine the ancestral allele in finite sites simulations.
#'
#' @param population The population that is marked as outgroup.
#' @export
#' @examples
#' # A simple finite sites model
#' model <- coal_model(c(4, 6, 1), 2, 10) +
#'    feat_outgroup(3) +
#'    feat_pop_merge(par_range('tau', 0.5, 2), 2, 1) +
#'    feat_pop_merge(par_expr('2*tau'), 3, 1) +
#'    feat_mutation(par_range('theta', 1, 10), model="HKY")
feat_outgroup <- function(population) {
  Feature_outgroup$new(population)
}


conv_to_ms_arg.Feature_outgroup <- function(feature, model) {
  stop("ms does not support outgroups")
}

conv_to_msms_arg.Feature_outgroup <- function(feature, model) {
  stop("msms does not support outgroups")
}

conv_to_seqgen_arg.Feature_outgroup <- ignore_par
