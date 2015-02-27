#' @importFrom R6 R6Class
sumstat_sfs_class <- R6Class('sumstat_SFS', inherit = sumstat,
  private = list(population = NULL),
  public = list(
    initialize = function(name, population, group=0) {
      assert_that(length(population) == 1)
      private$population <- population
      super$initialize(name, group)
    },
    calculate = function(seg_sites, files, model) {
      if ('all' %in% private$population) {
        individuals <- 1:sum(get_sample_size(model))
      } else {
        individuals <- get_population_indiviuals(model, private$population)
      }
      as.vector(calc_jsfs(seg_sites, individuals, numeric()))
    }
  )
)

#' Calculates the Site Frequency Spectrum from simulations
#'
#' @inheritParams sumstat_file
#' @param population Either the number of a population for with the SFS is
#'   calculated, or \code{'all'} to calculate the combined SFS of all
#'   populations.
#' @export
sumstat_sfs <- function(name='sfs', population='all', group = 0) {
  sumstat_sfs_class$new(name, population, group)
}
