#' @importFrom R6 R6Class
Simulator <- R6Class("Simulator",
  private = list(
    name = "TEMPLATE",
    features = NULL,
    sumstats = NULL,
    priority = 50
  ),
  public = list(
    simulate = function() stop("virtual method"),
    get_cmd = function() stop("virtual method"),
    get_name = function() private$name,
    get_features = function() private$features,
    get_sumstats = function() private$sumstats,
    get_priority = function() private$priority,
    initialize = function() {}
  )
)

is_simulator <- function(simulator) "Simulator" %in% class(simulator)


# Keep a user modifiable list of availible simulation programs in a private
# enviroment
if (!exists("simulators")) simulators <- new.env()

register_simulator <- function(simulator) {
  sim <- simulator$new()
  assert_that(is_simulator(sim))
  simulators[[sim$get_name()]] <- sim
}

get_simulator <- function(name) {
  simulators[[name]]
}
