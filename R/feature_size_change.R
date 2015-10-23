size_change_class <- R6Class("size_change", inherit = feature_class,
  public = list(
    print = function() {
      cat("Instanious size change to", print_par(private$rate),
          "* N0 in population", self$get_population(),
          "at time", print_par(self$get_time()), "\n")
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
#' then use \link{feat_growth}.
#'
#' @param new_size A \code{\link{parameter}} giving the new size of the
#'   population, as a factor of N0.
#' @param population The number of the population whichs size changes.
#'          Can also be set to "all". Then the size changes applies to all
#'          populations.
#' @param time The time at which the populations size is changed.
#' @return A feature which can be added to a model.
#' @export
#' @examples
#' # A model with one smaller population
#' model <- coal_model(c(20,37), 88) +
#'   feat_size_change(.1, 2, time="1")
feat_size_change <- function(new_size, population = 1, time = "0") {
  size_change_class$new(new_size, population, time)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.size_change <- function(feature, model) {
  all_pops <- feature$get_population() == "all" ||
    (feature$get_population() == 1 && length(get_populations(model)) == 1)
  present <- feature$get_time() == "par(0)"

  if (present) {
    if (all_pops) cmd <- "-N"
    else cmd <- "-n"
  } else {
    if (all_pops) cmd <- "-eN"
    else cmd <- "-en"
  }

  paste0(cmd, "', ",
         ifelse(present, "", paste(feature$get_time(), ", ")),
         ifelse(all_pops, "", paste(feature$get_population(), ", ")),
         feature$get_rate(), ", '")
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.size_change <- conv_to_ms_arg.size_change

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.size_change <- conv_to_ms_arg.size_change

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.size_change <- ignore_par
