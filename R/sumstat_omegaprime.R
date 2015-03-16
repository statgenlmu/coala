#' @importFrom R6 R6Class
SumstatOmegaPrime <- R6Class('SumstatOmegaPrime', inherit = Sumstat, #nolint
  private = list(
    population = NULL,
    req_segsites = TRUE
  ),
  public = list(
    initialize = function(name, population) {
      assert_that(is.numeric(population))
      assert_that(length(population) == 1)
      private$population <- population
      super$initialize(name)
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
sumstat_omegaprime  <- function(name='omegaprime', population=1) {
  SumstatOmegaPrime$new(name, population) #nolint
}
