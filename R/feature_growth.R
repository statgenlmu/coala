Feature_growth <- R6Class("Feature_growth", inherit = Feature,
  public = list(
    print = function() {
      cat("Exponential growth/decline with rate", self$get_par(),
          "in population", self$get_population(),
          "starting at time", self$get_time_point(), "\n")
    }
  )
)

#' Adds an exponential growth or decline of the size of one
#' population to a model.
#'
#' This function changes the growth factor of a population at given
#' point in time ('at.time'). This factor than applies to the time
#' interval farther into the past from this point.
#'
#' The population size changes by factor exp(-alpha*t), where alpha
#' is the growth parameter and t is the time since the growth has
#' started. Hence, for positive alpha, the population will 'decline
#' backwards in time' or grow forwards in time. Similar, will decline
#' in forwards time for a negative value of alpha.
#'
#' If you want to add an instantaneous change of the population size,
#' then use the \link{dm.addSizeChange} function.
#'
#' @param rate A \code{\link{Parameter}} stating the rate of the change.
#' @param population The population which starts to grow or decline.
#' @param at.time The time at which the population starts to grow or decline.
#' @return    The demographic model with a size change.
#' @export
#' @examples
#' dm <- coal_model(c(20,37), 88) +
#'   feat_growth(par_range('alpha', 0.1, 2), population=2, at.time="0")
feat_growth <- function(rate, population, at.time="0") {
  Feature_growth$new("growth", rate, pop_source=population, time_point=at.time)
}
