#' @importFrom R6 R6Class
stat_sfs_class <- R6Class("stat_sfs", inherit = sumstat_class,
  private = list(
    population = NULL,
    req_segsites = TRUE
  ),
  public = list(
    initialize = function(name, population) {
      assert_that(length(population) == 1)
      private$population <- population
      super$initialize(name)
    },
    calculate = function(seg_sites, trees, files, model) {
      if ("all" %in% private$population) {
        individuals <- 1:sum(get_sample_size(model))
      } else {
        individuals <- get_population_indiviuals(model, private$population)
      }
      sfs <- as.vector(calc_jsfs(seg_sites, individuals, numeric()))
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
sumstat_sfs <- function(name="sfs", population="all") {
  stat_sfs_class$new(name, population)
}
