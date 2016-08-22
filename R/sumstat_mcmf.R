#' @importFrom R6 R6Class
stat_mcmf_class <- R6Class("stat_mcmf", inherit = sumstat_class,
  private = list(
    population = NULL,
    req_segsites = TRUE,
    expand_mcmf = 1
  ),
  public = list(
    initialize = function(name, population, transformation, expand_mcmf) {
      assert_that(is.numeric(population))
      assert_that(length(population) == 1)
      private$population <- population
      private$expand_mcmf <- expand_mcmf
      super$initialize(name, transformation)
    },
    calculate = function(seg_sites, trees, files, model) {
      ploidy <- ifelse(is_unphased(model), get_ploidy(model), 1)
      #browser()
      calc_mcmf(seg_sites,
                get_population_individuals(model,
                                          private$population,
                                          haploids = (ploidy == 1)),
                get_locus_length_matrix(model),
                private$expand_mcmf,
                has_trios(model),
                ploidy)
    }
  )
)

#' Summary Statistic: Most Common Mutation's Frequency
#'
#' This summary statistic calculates the observed frequency
#' of the mutational pattern that is observed most often in
#' the data.
#'
#' The expand_mcmf = TRUE calculates only the mcmf per locus.
#' The expand_mcmf = 2 versions adds the frequency of
#' derived alleles in the most frequently observed mutational
#' pattern. The expand_mcmf = 3 also calculates the percentage
#' of positions that are polymorpic.
#'
#' @param name The name of the summary statistic. When simulating
#' a model, the value of the statistics are written to an entry
#' of the returned list with this name. Summary statistic names
#' must be unique in a model.
#' @param population The population for which the statistic is
#' calculated. Can also be "all" to calculate it from all populations.
#' @param transformation An optional function for transforming
#' the results of the statistic. If specified, the results of the
#' transformation are returned instead of the original values.
#' @param expand_mcmf The type of mcmf to be used. See Details
#' @return A numeric vector containing MCMF for each locus.
#'   \describe{
#'    \item{mcmf}{The observed frequency of the mutational pattern
#'      that is observed most often in the data.}
#'    \item{bal}{The frequency of derived alleles in the most
#'      frequently observed mutational pattern.}
#'    \item{perc_poly}{The percentage of positions that are polymorpic.}
#'   }
#'
#' @template summary_statistics
#' @examples
#' # Calculate MCMF for a panmictic population
#' model <- coal_model(10, 2) +
#'   feat_mutation(50) +
#'   sumstat_mcmf()
#' simulate(model)
#' @export
sumstat_mcmf  <- function(name = "mcmf", population = 1,
                          transformation = identity, expand_mcmf = TRUE) {
  stat_mcmf_class$new(name, population, transformation, expand_mcmf)
}
