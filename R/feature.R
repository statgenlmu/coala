self <- private <- super <- NULL # Mute warnings about R6 objects

Feature <- R6Class("Feature",
  private = list(
    parameter = list(),
    population = NULL,
    time = NULL,
    rate = NA,
    set_population = function(population, expected_length = 1) {
      assert_that(is.numeric(population))
      assert_that(length(population) == expected_length)
      private$population = population
    },
    add_parameter = function(parameter, required = TRUE, add_par = TRUE) {
      expr <- NA
      if (is.numeric(parameter) && length(parameter) == 1) {
        expr <- parameter
      } else if (is.character(parameter) && length(parameter) == 1) {
        expr <- parameter
      } else if (is.par(parameter)) {
        idx <- as.character(length(private$parameter) + 1)
        private$parameter[[idx]] <- parameter
        expr <- parameter$get_expression()
      } else if (is.null(parameter) || is.na(parameter)) {
        if (required) stop("Missing value for required parameter")
        return(NA)
      } else stop("Unexpected type of parameter")
      if (!add_par) return(expr)
      paste0("par(", expr, ")")
    }
  ),
  public = list(
    initialize = function(rate, population, time) {
      if (!missing(rate)) private$rate = private$add_parameter(rate)
      if (!missing(time)) private$time = private$add_parameter(time)
      if (!missing(population)) private$set_population(population)
    },
#     initialize = function(type, parameter,
#                           pop_source=NA, pop_sink=NA,
#                           time_point=NA,
#                           variance=0, zero_inflation=0) {
#
#       # Add the time point, which might also be a parameter
#       if (is.numeric(time_point)) time_point <- as.character(time_point)
#       else if (is.par(time_point)) {
#         self$add_parameter(time_point)
#         time_point <- as.character(time_point$get_expression())
#       }
#       else if (is.character(time_point) | is.na(time_point)) NULL
#       else stop("Unexpected type of argument 'time_point'")
#       private$time_point <- time_point
#
#       # Add the primary parameter (if any)
#       if (is.numeric(parameter)) par_expr <- as.character(parameter)
#       else if (is.character(parameter)) par_expr <- parameter
#       else if (is.par(parameter)) {
#         self$add_parameter(parameter)
#         par_expr <- as.character(parameter$get_expression())
#       }
#       else stop("Unexpected type of argument 'parameter'")
#
#       if (variance != 0) {
#         private$inter_locus_var <- TRUE
#         par_expr <- paste0('rgamma(1, ', par_expr, '^2/', variance,
#                             ', ', par_expr, '/', variance, ')')
#       }
#
#       if (zero_inflation != 0) {
#         private$inter_locus_var <- TRUE
#         par_expr <- paste0('ifelse(locus <= ',
#                             zero_inflation, ' * locus_number',
#                             ', 0, ', par_expr, ')')
#       }
#       private$par <- par_expr
#
#       private$population <- pop_source
#       private$feature_table <- create_feature_table(type, par_expr, pop_source,
#                                                     pop_sink, time_point)
#     },
    get_parameters = function() private$parameter,
    get_population = function() private$population,
    get_time = function() private$time,
    reset_parameters = function() private$parameter <- list(),
    print = function() {
      cat("Feature of type", private$feature_table$type[1], "\n")
    },
    get_call = function() private$call,
    get_rate = function() private$rate
  )
)

is.feature <- function(feature) {
  'Feature' %in% class(feature)
}

ignore_par <- function(feature, model) ""

print_par <- function(par) paste0("`", substr(par, 5, nchar(par)-1), "`")
