Feature_growth <- R6Class("Feature_growth", inherit = Feature,
  public = list(
    print = function() {
      cat("Exponential growth/decline with rate", print_par(private$rate),
          "in population", self$get_population(),
          "starting at time", print_par(self$get_time()), "\n")
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
#' then use the \link{model.addSizeChange} function.
#'
#' @param rate A \code{\link{Parameter}} stating the rate of the change.
#' @param population The population which starts to grow or decline.
#' @param time The time at which the population starts to grow or decline.
#' @return    The demographic model with a size change.
#' @export
#' @examples
#' model <- coal_model(c(20,37), 88) +
#'   feat_growth(par_range('alpha', 0.1, 2), population=2, time="0")
feat_growth <- function(rate, population, time="0") {
  Feature_growth$new(rate, population, time)
}


conv_to_ms_arg.Feature_growth <- function(feature, model) {
  all_pops <- feature$get_population() == "all" ||
    (feature$get_population() == 1 && length(get_populations(model)) == 1)
  present <- feature$get_time() == "par(0)"

  if (present) {
    if (all_pops) cmd <- "-G"
    else cmd <- "-g"
  } else {
    if (all_pops) cmd <- "-eG"
    else cmd <- "-eg"
  }

  paste0(cmd, "', ",
         ifelse(present, "", paste(feature$get_time(), ", ")),
         ifelse(all_pops, "", paste(feature$get_population(), ", ")),
         feature$get_rate(), ", '")
}

conv_to_msms_arg.Feature_growth <- conv_to_ms_arg.Feature_growth
conv_to_scrm_arg.Feature_growth <- conv_to_ms_arg.Feature_growth
conv_to_seqgen_arg.Feature_growth <- ignore_par
