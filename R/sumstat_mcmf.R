#' @importFrom R6 R6Class
stat_mcmf_class <- R6Class("stat_mcmf", inherit = sumstat_class,
  private = list(
    population = NULL,
    req_segsites = TRUE
  ),
  public = list(
    initialize = function(name, population, transformation) {
      assert_that(is.numeric(population))
      assert_that(length(population) == 1)
      private$population <- population
      super$initialize(name, transformation)
    },
    calculate = function(seg_sites, trees, files, model) {
      ploidy <- ifelse(is_unphased(model), get_ploidy(model), 1)
      calc_mcmf(seg_sites,
                get_population_indiviuals(model,
                                          private$population,
                                          haploids = (ploidy == 1)),
                has_trios(model),
                ploidy)
    }
  )
)

#' The MCMF Summary Statistic
#'
#' When adding this to a model, the "Most Common Mutation Frequecy" summary
#' statistic is calculated from the simulation output.
#'
#' @inheritParams sumstat_four_gamete
#' @export
sumstat_mcmf  <- function(name = "mcmf", population = 1,
                          transformation = identity) {
  stat_mcmf_class$new(name, population, transformation)
}
