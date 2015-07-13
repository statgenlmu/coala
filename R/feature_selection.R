selection_class <- R6Class("selection", inherit = feature_class,
  private = list(strength_AA = NA, strength_Aa = NA, additive = FALSE),
  public = list(
    initialize = function(strength_AA, strength_Aa, population, time,
                          additive = FALSE) {
      if (additive) {
        private$strength_AA <- private$add_parameter(strength_AA)
        private$additive <- TRUE
      } else {
        private$strength_AA <- private$add_parameter(strength_AA)
        private$strength_Aa <- private$add_parameter(strength_Aa)
      }

      private$set_population(population)
      private$time <- private$add_parameter(time)
    },
    get_strength_AA = function() private$strength_AA,
    get_strength_Aa = function() private$strength_Aa,
    is_additive = function() private$additive,
    print = function() {
      if (private$additive) {
        cat("Additive selection with strength ", print_par(private$strength_AA))
      } else {
        cat("Selection with strength ", print_par(private$strength_AA),
            "(AA) and ", print_par(private$strength_Aa), "(Aa)")
      }
      cat(" in population", self$get_population(),
          "starting at time", print_par(self$get_time()), "\n")
    }
  )
)

#' Adds positive selection to a model
#'
#' @param population The population in which the allele is selected.
#' @param time The time at which the selection starts.
#' @param strength_AA The selection strength for the selected homozygote
#' @param strength_Aa The selection strength for the heterozygote.
#' @param strength_A This sets the strength for the selected allele in an
#'   haploid model. \code{strength_AA} and \code{strength_Aa} are ignored
#'   when this is argument is given.
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
feat_selection <- function(strength_AA, strength_Aa, strength_A,
                           population = 1, time) {
  if (!missing(strength_A)) {
    return(selection_class$new(strength_A,
                               population = population,
                               time = time,
                               additive = TRUE))
  }
  selection_class$new(strength_AA, strength_Aa, population, time)
}


#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.selection <- function(feature, model) {
  stop("selection is not supported", call. = FALSE)
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.selection <- conv_to_ms_arg.selection

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.selection <- function(feature, model) {
  n_pop <- length(get_populations(model))
  start_freq <- rep(0, n_pop)
  start_freq[feature$get_population()] <- 0.0005
  if (feature$is_additive()) {
    strength <- paste0("-SA',", feature$get_strength_AA(), ", '")
  } else {
    strength <- paste0("-SAA',", feature$get_strength_AA(), ", '",
                       "-SAa',", feature$get_strength_Aa(), ", '")
  }
  paste0("-SI', ", feature$get_time(), ", '",
         n_pop, " " , paste(start_freq, collapse = " "), " ",
         strength,
         "-Sp 0.5 -SForceKeep -N 10000 ")
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.selection <- ignore_par
