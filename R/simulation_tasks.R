#' @importFrom digest digest
simulation_task_class <- R6Class("simulation_task",
  private = list(
    n_loci = 0,
    args = list(),
    simulator = NULL,
    args_hash = ""
  ),
  active = list(
    locus_number = function(n) {
      if (missing(n)) return(private$n_loci)
      assert_that(is.number(n))
      names(n) <- NULL
      private$n_loci <- n
    }
  ),
  public = list(
    initialize = function(simulator, locus_number, ...) {
      assert_that(is_simulator(simulator))
      self$locus_number <- locus_number
      private$args <- list(...)
      # Repect simulator when hashing:
      private$args[["_simulator"]] <- simulator$get_name()
      private$args_hash <- digest(private$args)
      private$simulator <- simulator
    },
    hash = function() private$args_hash,
    get_arg = function(arg) {
      if (!any(arg == names(private$args))) stop("Argument ", arg, " not found")
      private$args[[arg]]
    },
    get_simulator = function() private$simulator,
    print = function() {
      cat("Simulation Task\n")
      cat("---------------\n")
      cat("Simulator: ", private$simulator$get_name(), "\n")
      cat("Locus Number: ", self$locus_number, "\n")
      for (i in seq_along(private$args)) {
        if (names(private$args)[i] == "_simulator") next
        cat(names(private$args)[i], ": ")
        print(private$args[[i]])
      }
      cat("Hash: ", self$hash(), "\n")
      cat("---------------\n")
    }
  )
)


is_simulation_task <- function(x) inherits(x, "simulation_task")

create_sim_task <- function(simulator, locus_number, ...) {
  simulation_task_class$new(simulator, locus_number, ...)
}


generate_sim_tasks <- function(model, pars) {
  simulator <- select_simprog(model)
  if (is.null(simulator)) stop("No simulator found")

  locus_groups <- seq_len(get_locus_group_number(model))
  tasks <- lapply(locus_groups, function(locus_group) {
    group_model <- create_group_model(model, locus_group)
    total_locus_number <- get_locus_number(group_model, group = 1, TRUE)
    if (has_variation(group_model)) {
      locus_number <- rep(1, total_locus_number)
    } else {
      locus_number <- total_locus_number
    }

    sim_tasks <- lapply(seq_along(locus_number), function(locus_id) {
      simulator$create_task(group_model, pars,
                            locus_number[locus_id],
                            locus_id)
    })

    reduce_sim_tasks(sim_tasks)
  })

  unlist(tasks, recursive = FALSE)
}


fill_cmd_template <- function(template, model, parameters,
                              locus_id, eval_pars = TRUE) {

  locus_length <- get_locus_length(model, group = 1)
  locus_number <- get_locus_number(model, 1, TRUE)

  tmp_env <- create_par_env(model, parameters,
                            locus_length = locus_length,
                            locus_id = locus_id,
                            locus_number = locus_number,
                            for_cmd = !eval_pars)

  paste(eval(parse(text = template), envir = tmp_env), collapse = " ")
}


reduce_sim_tasks <- function(tasks) {
  hashes <- vapply(tasks, function(x) x$hash(), character(1))
  locus_numbers <- vapply(tasks, function(x) x$locus_number, numeric(1))

  locus_numbers_reduced <- by(locus_numbers, hashes, sum)
  lapply(tasks[!duplicated(hashes)], function(x) {
    x$locus_number <- locus_numbers_reduced[x$hash()]
    x
  })
}


combine_results <- function(sim_results) {
  results <- list()

  if (!is.null(sim_results[[1]]$seg_sites)) {
    results$seg_sites <-
      do.call(c, lapply(sim_results, function(x) x$seg_sites))
  }

  if (!is.null(sim_results[[1]]$trees)) {
    results$trees <- do.call(c, lapply(sim_results, function(x) x$trees))
  }

  if (!is.null(sim_results[[1]]$files)) {
    results$files <- do.call(c, lapply(sim_results, function(x) x$files))
  }

  results$cmds <- do.call(list, lapply(sim_results, function(x) x$cmd))
  results$simulators <- do.call(list, lapply(sim_results,
                                             function(x) x$simulator))
  results
}
