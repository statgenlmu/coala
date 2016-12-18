#' @importFrom R6 R6Class
stat_mcmf_class <- R6Class("stat_mcmf", inherit = sumstat_class,
  private = list(
    population = NULL,
    req_segsites = TRUE,
    expand_mcmf = FALSE,
    type_expand = 1
  ),
  public = list(
    initialize = function(name, population, transformation,
                          expand_mcmf, type_expand) {

      assert_that(identical(population, "all") || is.numeric(population))

      assert_that(length(population) == 1)
      assert_that(is.logical(expand_mcmf))
      assert_that(is.numeric(type_expand))
      private$population <- population
      private$expand_mcmf <- expand_mcmf
      private$type_expand <- type_expand
      super$initialize(name, transformation)
    },
    calculate = function(seg_sites, trees, files, model, sim_tasks = NULL) {
      ploidy <- get_samples_per_ind(model)
      individuals <- get_population_individuals(model,
                                                private$population,
                                                haploids = (ploidy == 1))
      locus_length <- get_locus_length_matrix(model)
      if (is.null(locus_length)) locus_length <- matrix(0, 0, 6)

      mcmf <- calc_mcmf(seg_sites = seg_sites,
                        individuals = individuals,
                        expand_mcmf = private$expand_mcmf,
                        type_expand = private$type_expand,
                        has_trios = has_trios(model),
                        ploidy = ploidy,
                        locus_length = locus_length)

      if (private$expand_mcmf == FALSE) {
        return(mcmf[, 1])
      }

      mcmf
    }
  )
)

#' Summary Statistic: Most Common Mutation's Frequency
#'
#' This summary statistic calculates the observed frequency
#' of the mutational pattern that is observed most often in
#' the data.
#'
#' The expand_mcmf = FALSE calculates the mcmf per locus
#' and returns a vector. The expand_mcmf = TRUE and type_expand = 1
#' returns the same results as the first column of a Matrix. The
#' expand_mcmf = TRUE and type_expand = 2 adds the frequency of
#' derived alleles in the most frequently observed mutational pattern
#' as a second column. The expand_mcmf = TRUE and type_expand = 3 adds
#' the percentage of positions that are polymorpic. When
#' expanded_mcmf = TRUE results are returned as a matrix.
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
#' @param expand_mcmf Whether to use or not the expanded MCMF. See Details
#' @param type_expand Specifies the type of expanded MCMF to be used. See Details
#' @return A numeric vector or matrix containing MCMF for each locus.
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
sumstat_mcmf  <- function(name = "mcmf",
                          population = 1,
                          transformation = identity,
                          expand_mcmf = FALSE,
                          type_expand = 1) {

  stat_mcmf_class$new(name,
                      population,
                      transformation,
                      expand_mcmf,
                      type_expand)
}
