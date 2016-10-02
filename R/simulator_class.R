#' @importFrom R6 R6Class
simulator_class <- R6Class("simulator",
  private = list(
    name = "TEMPLATE",
    priority = 50
  ),
  public = list(
    # Methods that simulators must implement
    create_task = function(model, pars, locus_number,
                           locus_id = 1,
                           eval_pars = TRUE) stop("virtual method"),
    simulate = function(model, sim_task) stop("virtual method"),

    # Infrastructure for simulators
    get_name = function() private$name,
    get_priority = function() private$priority,
    get_info = function() c(name = private$name),
    initialize = function(priority) {
      assert_that(is.number(priority))
      private$priority <- priority
    }
  )
)


is_simulator <- function(simulator) inherits(simulator, "simulator")


# Keep a user modifiable list of available simulation programs in a private
# environment
if (!exists("simulators")) simulators <- new.env()

register_simulator <- function(simulator) {
  assert_that(is_simulator(simulator))
  simulators[[simulator$get_name()]] <- simulator
}

get_simulator <- function(name) {
  simulators[[name]]
}





#' Returns the available simulators
#'
#' This functions returns the usable simulators
#'
#' @export
#' @examples
#' list_simulators()
list_simulators <- function() {
  simulators <- do.call(rbind, lapply(ls(simulators), function(simulator) {
    info <- get_simulator(simulator)$get_info()
    name <- info[["name"]]
    info <- info[-1]
    pars <- paste(names(info), ":", info, collapse = ", ")
    data.frame(name = name,
               priority = get_simulator(simulator)$get_priority(),
               info = pars)
  }))
  simulators[order(simulators$priority, decreasing = TRUE), ]
}


test_simulator_class <- R6Class("test", inherit = simulator_class,
  private = list(
    name = "test",
    priority = -99999
  ),
  public = list(
    initialize = function() {}, #nolint
    create_cmd_template = function(model) {
      pars <- paste0("par(", get_par_names(model), ")")
      paste0("c('sum(', ", paste(pars, collapse = ", '+', "), ", ')')")
    },
    get_cmd = function(model) {
      template <- self$create_cmd_template(model)
      cmd <- fill_cmd_template(template, model, NULL, 1, eval_pars = FALSE)
      paste("test",
            cmd[1, "locus_number"],
            cmd[1, "command"])
    },
    simulate = function(model, cmd) {
      eval()
    }
  )
)

test_simulator <- function() test_simulator_class$new()

test_feature_class <- R6Class("test_feature", inherit = feature_class)

test_model <- function(locus_number = 10) {
  coal_model(5, locus_number) +
    test_feature_class$new() +
    par_named("a") +
    par_named("b")
}
