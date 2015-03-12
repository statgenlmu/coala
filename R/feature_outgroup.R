Feature_outgroup <- R6Class("Feature_outgroup", inherit = Feature,
  public = list(
    print = function() {
      cat("Outgroup: Population", self$get_par(), "\n")
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
#' dm <- coal_model(c(4, 6, 1), 2, 10) +
#'    feat_outgroup(3) +
#'    feat_pop_merge(par_range('tau', 0.5, 2), 2, 1) +
#'    feat_pop_merge(par_expr('2*tau'), 3, 1) +
#'    feat_mutation(par_range('theta', 1, 10), model="HKY")
feat_outgroup <- function(population) {
  Feature_outgroup$new("outgroup", par_const(population))
}
