Feature_size_change <- R6Class("Feature_size_change", inherit = Feature,
  private = list(rate = NA),
  public = list(
    initialize = function(new_size, population, time) {
      private$rate = private$add_parameter(new_size)
      private$time = private$add_parameter(time)
      private$set_population(population)
    },
    get_rate = function() private$rate,
      print = function() {
        cat("Instanious size change to", private$rate,
            "* N0 in population", self$get_population(),
            "at time", self$get_time(), "\n")
    }
  )
)


#' Adds an instantaneous change of the population size of one
#' population to a model.
#'
#' This function changes the effective population size of one
#' population. The change is performed at a given time point
#' and applies to the time interval further on into
#' the past from this point. The population size is set to a
#' fraction of N0, the present day size of population one.
#'
#' If you want to add a slow, continuous change over some time,
#' then use the \link{model.addGrowth} function.
#'
#' @param new_size A \code{\link{parameter}} giving the new size of the
#'   population, as a factor of N0.
#' @param population The number of the population whichs size is changed.
#' @param time The time point at which the size changes.
#' @return A feature which can be added to a model.
#' @export
#' @examples
#' # A model with one smaller population
#' model <- coal_model(c(20,37), 88) +
#'   feat_size_change(.1, 2, time="1")
feat_size_change <- function(new_size, population, time="0") {
  Feature_size_change$new(new_size, population, time)
}

conv_to_ms_arg.Feature_size_change <- function(feature, model) {
  paste0("-en\", ", feature$get_time(), ", \"",
         feature$get_population(), "\", ",
         feature$get_rate(), ", \"")
}

conv_to_msms_arg.Feature_size_change <- ignore_par
conv_to_seqgen_arg.Feature_size_change <- ignore_par
