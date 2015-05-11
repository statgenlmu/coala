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
