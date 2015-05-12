#' @importFrom R6 R6Class
Simulator <- R6Class("Simulator",
  private = list(
    name = "TEMPLATE",
    priority = 50
  ),
  public = list(
    simulate = function() stop("virtual method"),
    get_cmd = function() stop("virtual method"),
    get_name = function() private$name,
    get_features = function() private$features,
    get_sumstats = function() private$sumstats,
    get_priority = function() private$priority,
    initialize = function() NULL
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


fill_cmd_template <- function(template, model, parameters,
                              locus_group, eval_pars = TRUE) {
  locus_length <- get_locus_length(model, group = locus_group)

  if (has_variation(model)) {
    locus_number <- rep(1, get_locus_number(model, locus_group, TRUE))
  } else {
    locus_number <- get_locus_number(model, locus_group, TRUE)
  }

  args <- sapply(locus_number, function(ln) {
    tmp_env <- create_par_env(model, parameters,
                              locus_length = locus_length,
                              for_cmd = !eval_pars)

    paste(eval(parse(text = template), envir = tmp_env), collapse = " ")
  })

  sim_cmds <- data.frame(locus_number = locus_number,
                         command = args,
                         stringsAsFactors = FALSE)

  reduce_sim_commands(sim_cmds)
}


reduce_sim_commands <- function(sim_commands) {
  if (nrow(sim_commands) == 1) return(sim_commands)
  grouped_commands <- unique(sim_commands[ , 2])
  if (length(grouped_commands) == nrow(sim_commands)) return(sim_commands)
  grouped_locus_number <- vapply(grouped_commands, function(cmd) {
    sum(sim_commands[sim_commands[ , 2] == cmd, 1])
  }, numeric(1))
  data.frame(locus_number = grouped_locus_number,
             command = grouped_commands,
             stringsAsFactors = FALSE)
}
