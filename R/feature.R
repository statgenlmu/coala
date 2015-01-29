Base_Object <- R6Class("CSR_OBJ")

Feature <- R6Class("Feature", inherit = Base_Object,
  private = list(
    feature_table = NULL,
    parameter = list(),
    group = 0,
    inter_locus_var = FALSE
  ),
  public = list(
    initialize = function(type, parameter,
                          pop_source=NA, pop_sink=NA,
                          time_point=NA, group=0,
                          variance=0, zero_inflation=0) {

      # Add the parameter
      if (is.numeric(parameter)) par_expr <- as.character(parameter)
      else if (is.character(parameter)) par_expr <- parameter
      else if (is.par(parameter)) {
        self$add_parameter(parameter)
        par_expr <- as.character(parameter$get_expression())
      }
      else stop("Unexpected type of argument 'parameter'")

      # Add the time point, which might also be a parameter
      if (is.numeric(time_point)) time_point <- as.character(time_point)
      else if (is.par(time_point)) {
        self$add_parameter(time_point)
        time_point <- as.character(time_point$get_expression())
      }
      else if (is.character(time_point) | is.na(time_point)) {}
      else stop("Unexpected type of argument 'time_point'")

      private$group <- group

      if (variance != 0) {
        private$inter_locus_var <- TRUE
        par_expr <- paste0('rgamma(1, ', par_expr, '^2/', variance,
                            ', ', par_expr, '/', variance, ')')
      }

      if (zero_inflation != 0) {
        private$inter_locus_var <- TRUE
        par_expr <- paste0('ifelse(locus <= ',
                            zero_inflation, ' * locus_number',
                            ', 0, ', par_expr, ')')
      }

      private$feature_table <- createFeatureTable(type, par_expr, pop_source,
                                                  pop_sink, time_point, group)
    },
    get_table = function() private$feature_table,
    get_parameters = function() private$parameter,
    get_group = function() private$group,
    get_inter_locus_var = function() private$inter_locus_var,
    add_parameter = function(parameter) {
      stopifnot(is.par(parameter))
      idx <- as.character(length(private$parameter) + 1)
      private$parameter[[idx]] <- parameter
    },
    add_feature = function(feature) {
      stopifnot(is.feature(feature))
      private$feature_table <- rbind(private$feature_table, feature$get_table())
      for (par in feature$get_parameters()) {
        self$add_parameter(par)
      }
    }
  )
)

is.feature <- function(feature) {
  'Feature' %in% class(feature)
}

