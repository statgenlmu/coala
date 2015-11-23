#' @importFrom R6 R6Class
stat_sfs_class <- R6Class("stat_sfs", inherit = sumstat_class,
  private = list(
    population = NULL,
    req_segsites = TRUE
  ),
  public = list(
    initialize = function(name, population, transformation) {
      assert_that(length(population) == 1)
      private$population <- population
      super$initialize(name, transformation)
    },
    calculate = function(seg_sites, trees, files, model) {
      individuals <- get_population_indiviuals(model, private$population)
      sfs <- as.vector(calc_jsfs(seg_sites, list(individuals)))
      sfs[c(-1, -length(sfs))]
    }
  )
)

#' Calculates the Site Frequency Spectrum from simulations
#'
#' @inheritParams sumstat_four_gamete
#' @param population Either the number of a population for with the SFS is
#'   calculated, or \code{"all"} to calculate the combined SFS of all
#'   populations.
#' @export
#' @examples
#' set.seed(50)
#' model <- coal_model(5, 2) + feat_mutation(5) + sumstat_sfs()
#' stats <- simulate(model)
#' print(stats$sfs)
sumstat_sfs <- function(name = "sfs", population = "all",
                        transformation = identity) {
  stat_sfs_class$new(name, population, transformation)
}
