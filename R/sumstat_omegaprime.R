#' @importFrom R6 R6Class
SumStat_OmegaPrime <- R6Class('SumStat_OmegaPrime', inherit = SumStat,
  private = list(population = NULL),
  public = list(
    initialize = function(name, population, group=0) {
      assert_that(is.numeric(population))
      assert_that(length(population) == 1)
      private$population <- population
      super$initialize(name, group)
    },
    calculate = function(seg_sites, files, model) {
      calc_omegaprime(seg_sites,
                      get_population_indiviuals(model, private$population))
    }
  )
)

#' Calculates the (experimental) Omega' Statistic
#'
#' @inheritParams sumstat_file
#' @export
sumstat_omegaprime  <- function(name='omegaprime', population=1, group = 0) {
  SumStat_OmegaPrime$new(name, population, group)
}
