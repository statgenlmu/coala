Feature_selection <- R6Class("Feature_selection", inherit = Feature,
  private = list(strength_AA = NA, strength_Aa = NA),
  public = list(
    initialize = function(strength_AA, strength_Aa, population, time) {
      private$strength_AA = private$add_parameter(strength_AA)
      private$strength_Aa = private$add_parameter(strength_Aa)
      private$set_population(population)
      private$time = private$add_parameter(time)
    },
    get_strength_AA = function() private$strength_AA,
    get_strength_Aa = function() private$strength_Aa,
    print = function() {
    cat("Selection with strength ", print_par(private$strength_AA),
        "(AA) and ", print_par(private$strength_Aa), "(Aa)",
        "in population", self$get_population(),
        "starting at time", print_par(self$get_time()), "\n")
    }
  )
)

#' Adds positiv selection to a model
#'
#' @param population The populaton in which the allele is selected.
#' @param time The time at which the selection starts.
#' @export
#' @examples
#' # Positive selection in population 2:
#' model <- coal_model(c(10, 13), 100) +
#'   feat_pop_merge(par_range('tau', .1, 2), 2, 1) +
#'   feat_selection(strength_AA=par_expr(2*s),
#'                  strength_Aa=par_range('s', 100, 2000),
#'                  population = 2,
#'                  time=par_expr(tau))
#'
feat_selection <- function(strength_AA, strength_Aa, population = 1, time) {
  Feature_selection$new(strength_AA, strength_Aa, population, time)
}

conv_to_ms_arg.Feature_selection <- function(feature, model) {
  stop("ms does not support selection")
}

conv_to_msms_arg.Feature_selection <- function(feature, model) {
  n_pop <- length(get_populations(model))
  start_freq <- rep(0, n_pop)
  start_freq[feature$get_population()] <- 0.0005
  paste0("-SI', ", feature$get_time(), ", '",
         n_pop, " " , paste(start_freq, collapse = " "), " ",
         "-SAA',", feature$get_strength_AA(), ", '",
         "-SAa',", feature$get_strength_Aa(), ", '",
         "-Sp 0.5 -SForceKeep -N 10000 ")
}

conv_to_seqgen_arg.Feature_selection <- ignore_par
