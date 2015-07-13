migration_class <- R6Class("migration", inherit = feature_class,
  private = list(rate = NA),
  public = list(
    initialize = function(rate, pop_from, pop_to, time, symmetric = FALSE) {
      super$initialize(rate = rate, time = time)

      if (symmetric) {
        private$population <- "all"
      } else {
        private$set_population(c(from = pop_from, to = pop_to), 2)
      }
    },
    print = function() {
      if (all(private$population == "all")) {
        cat("Symmetric migration")
      } else {
        cat("Migration from pop", private$population,
            "to pop", private$pop_to)
      }
      cat(" with rate", print_par(private$rate),
          "starting at time", print_par(self$get_time()), "\n")
    }
  )
)


#' Add migration/gene flow between two populations to a demographic model
#'
#' This function adds the assumption to the model that some individuals
#' 'migrate' from one sub-population to another, i.e. they leave the one
#' and become a member of the other. This is usually used to model ongoing
#' gene flow through hybridisation after the populations separated.
#'
#' You can enter a time (\code{time}) at which the migration is
#' assumed to start (looking backwards in time). From that time on, a
#' fixed number of migrants move from population \code{pop_from} to
#' population \code{pop_to} each generation. This number is given via this
#' feature's parameter, which equals 4*N0*m,  where m is the
#' fraction of \code{pop_to}that is replaced with migrants each generation.
#' If \code{pop_to} has also size Ne, than this is just the
#' expected number of individuals that migrate each generation.
#'
#' You can add different mutation rates at different times to your model.
#' Then each rate will be used for the period from its time point to
#' the next. Migration from and to an population always ends with the
#' speciation event in which the population is created.
#'
#' @param rate  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "M" will use
#'            an parameter with name M that you have previously
#'            created. You can also use R expression here, so "2*M"
#'            or "5*M+2*tau" (if tau is another parameter) will also
#'            work (also this does not make much sense).
#' @param pop_from The population from which the individuals leave.
#' @param pop_to The population to which the individuals move.
#' @param symmetric Use the rate between all pairs of populations.
#' @param time The time point at which the migration with this rate starts.
#' @export
#'
#' @examples
#' # Asymmetric migration for two populations
#' model <- coal_model(c(25,25), 100) +
#'   feat_migration(0.5, 1, 2) +
#'   feat_migration(0.75, 2, 1)
#'
#' # Symmetric Migration
#' model <- coal_model(c(25,25), 100) +
#'   feat_migration(par_range('m', 0.1, 2), symmetric=TRUE)
feat_migration <- function(rate, pop_from = NULL, pop_to = NULL,
                           symmetric = FALSE, time = "0") {
  if (symmetric) {
    return(migration_class$new(rate, time = time, symmetric = TRUE))
  } else {
    return(migration_class$new(rate, pop_from, pop_to, time))
  }
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.migration <- function(feature, model) {
  if (all(feature$get_population() == "all")) {
    return( paste0("-eM', ", feature$get_time(), ", ",
                    feature$get_rate(), ", '"))
  }
  paste0("-em', ", feature$get_time(), ", ",
         feature$get_population()[1], ", ",
         feature$get_population()[2], ", ",
         feature$get_rate(), ", '")
}

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.migration <- conv_to_ms_arg.migration

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.migration <- conv_to_ms_arg.migration

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.migration <- ignore_par
