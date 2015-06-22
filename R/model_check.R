#' Check which simulator can simulate a model
#'
#' This function allows to check which of the available simulators can
#' simulate a given model. It also states the problems for the ones that
#' are incompatible with the model.
#'
#' @param model The model that is checked
#' @export
check_model <- function(model) {
  force(model)
  for (simprog_name in ls(simulators)) {
    cat(simprog_name, ": ")
    simprog <- get_simulator(simprog_name)

    cmd <- try(simprog$get_cmd(model), silent = TRUE)

    if (inherits(cmd, "try-error")) {
      cat(cmd, "\n")
    } else {
      cat("OK\n\n")
    }
  }
}
