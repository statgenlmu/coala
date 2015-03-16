#' @importFrom R6 R6Class
SumstatJsfs <- R6Class('SumstatJsfs', inherit = Sumstat, #nolint
  private = list(
    populations = NULL,
    req_segsites = TRUE
  ),
  public = list(
    initialize = function(name, populations) {
      assert_that(is.numeric(populations))
      assert_that(length(populations) == 2)
      private$populations <- populations
      super$initialize(name)
    },
    calculate = function(seg_sites, files, model) {
      calc_jsfs(seg_sites,
                get_population_indiviuals(model, private$populations[1]),
                get_population_indiviuals(model, private$populations[2]))
    }
  )
)

#' Calculates the Joint Site Frequency Spectrum from simulations
#'
#' @param name The name of the summary statistic. When simulating a model,
#'   the value of the statistics are written to an entry of the returned list
#'   with this name. Summary statistic names must be unique in a model.
#' @param population The population for which the statistic is calculated.
#' @export
sumstat_jsfs <- function(name='jsfs', populations=c(1,2)) {
  SumstatJsfs$new(name, populations) #nolint
}
