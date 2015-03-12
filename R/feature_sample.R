Feature_sample <- R6Class("Feature_sample", inherit = Feature,
  public = list(
    print = function() {
      cat("Sampling of", self$get_par(),
          ifelse(self$get_par() == 1, "haploid", "haploids"),
          "in population", self$get_population(),
          #"at time", self$get_time_point(),
          "\n")
    }
  )
)

#' Creates a feature that represents the sampling from one population
#'
#' @param size The number of individuals that are sampled.
#' @param population The population from with the indidivuals are sampled
#' @param time_point The time at which the sample is taken.
#' @return The feature, which can be added to a model using `+`.
#' @export
feat_sample <- function(size, population, time_point='0') {
  assert_that(is.numeric(size))
  if (time_point != '0')
    stop("Samples at time different from 0 at not supported at the moment")
  Feature_sample$new('sample',
                     size,
                     pop_source=population,
                     time_point=time_point)
}
