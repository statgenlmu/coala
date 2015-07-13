#' @importFrom R6 R6Class
stat_omega_prime_class <- R6Class("stat_omega_prime", inherit = sumstat_class,
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
    calculate = function(seg_sites, trees, files, model) {
      calc_omegaprime(seg_sites,
                      get_population_indiviuals(model, private$population))
    }
  )
)

#' Calculates the (experimental) Omega" Statistic
#'
#' @inheritParams sumstat_four_gamete
sumstat_omegaprime  <- function(name="omegaprime", population=1) {
  stat_omega_prime_class$new(name, population) #nolint
}
