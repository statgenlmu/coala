Feature_selection <- R6Class("Feature_selection", inherit = Feature,
  private = list(strength_AA = NA, strength_Aa = NA),
  public = list(
    initialize = function(strength_AA, strength_Aa, population, time) {
      private$strength_AA = private$add_parameter(strength_AA)
      private$strength_Aa = private$add_parameter(strength_Aa)
      private$set_population(population)
      private$time = private$add_parameter(time)
    },
    print = function() {
    cat("Exponential growth/decline with rate", private$rate,
      "in population", self$get_population(),
      "starting at time", self$get_time(), "\n")
    }
  )
)

#' Adds positiv selection to a model
#'
#' @param population The populaton in which the allele is selected.
#' @param at_time The time at which the selection starts.
#' @export
#' @examples
#' # Positive selection in population 2:
#' model <- coal_model(c(10, 13), 100) +
#'   feat_pop_merge(par_range('tau', .1, 2), 2, 1) +
#'   feat_selection(strength_AA=par_expr(2*s), strength_Aa=par_range('s', 100, 2000),
#'                  population = 2, at_time=par_expr(tau))
#'
feat_selection <- function(strength_AA, strength_Aa, population, at_time) {
  Feature_selection$new(strength_AA, strength_Aa, population, at_time)
}

conv_to_ms_arg.Feature_selection <- function(feature, model) {
  stop("ms does not support selection")
}

conv_to_msms_arg.Feature_growth <- function(feature, model) {
  #     if (type == "selection") {
  #       cmd <- c(cmd, '"-SI"', ',', feat['time.point'], ',',
  #                length(get_sample_size(model, for_sim = TRUE)), ',')
  #       start_freq <- rep(0, length(get_sample_size(model, for_sim = TRUE)))
  #       start_freq[ as.integer(feat['pop.source']) ] <- 0.0005
  #       cmd <- c(cmd, paste0('"', paste(start_freq, collapse=' '), '"'), ',')
  #
  #       cmd <- c(cmd, '"-N 10000"', ',')
  #
  #       s_AA <- search_feature(model, 'selection_AA',
  #                              pop.source = feat['pop.source'],
  #                              time.point = feat['time.point'])$parameter
  #       stopifnot(length(s_AA) == 1)
  #       cmd <- c(cmd, '"-SAA"', ',', s=s_AA, ',')
  #
  #       s_Aa <- search_feature(model, 'selection_Aa',
  #                             pop.source = feat['pop.source'],
  #                             time.point = feat['time.point'])$parameter
  #       stopifnot(length(s_Aa) == 1)
  #       cmd <- c(cmd, '"-SAa"', ',', s=s_Aa, ',')
  #       cmd <- c(cmd, '"-Sp 0.5"', ',', '"-SForceKeep"', ',')
  #       cmd <- c(cmd, '"-threads 1"', ',')
  #     }
  #   }

}

conv_to_seqgen_arg.Feature_growth <- ignore_par
