SumStat <- R6Class('SumStat', inherit = Base_Object,
  private = list(name='', population=0, group=0),
  public = list(
    initialize = function(name, population=0, group=0) {
      private$name <- name
      private$population <- population
      private$group <- group
    },
    get_name = function() private$name,
    get_population = function() private$population,
    get_group = function() private$group
  )
)


is.sum_stat <- function(sum_stat) 'SumStat' %in% class(sum_stat)
