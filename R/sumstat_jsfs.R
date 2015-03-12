#' @importFrom R6 R6Class
sumstat_jsfs_class <- R6Class('sumstat_Jsfs', inherit = sumstat,
  private = list(populations = NULL),
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
#' @inheritParams sumstat_file
#' @export
sumstat_jsfs <- function(name='jsfs', populations=c(1,2)) {
  sumstat_jsfs_class$new(name, populations)
}
