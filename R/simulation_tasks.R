generate_sim_tasks <- function(model, parameters) {
  simulator <- select_simprog(model)
  if (is.null(simulator)) stop("No simulator found")

  locus_groups <- seq_len(get_locus_group_number(model))
  tasks <- do.call(rbind, lapply(locus_groups, function(i) {
      cmd_template <- simulator$create_cmd_template(model)
      cmds <- fill_cmd_template(cmd_template, model, parameters, i)
      data.frame(simulator = simulator$get_name(), cmds,
                 stringsAsFactors = FALSE)
  }))

  tasks
}


fill_cmd_template <- function(template, model, parameters,
                              locus_group, eval_pars = TRUE) {

  locus_length <- get_locus_length(model, group = locus_group)
  total_locus_number <- get_locus_number(model, locus_group, TRUE)

  if (has_variation(model)) {
    locus_number <- rep(1, total_locus_number)
  } else {
    locus_number <- total_locus_number
  }
  locus_id <- seq(along = locus_number)

  args <- vapply(locus_id, function(l_id) {
    tmp_env <- create_par_env(model, parameters,
                              locus_length = locus_length,
                              locus_id = l_id,
                              locus_number = total_locus_number,
                              for_cmd = !eval_pars)

    paste(eval(parse(text = template), envir = tmp_env), collapse = " ")
  }, character(1L))

  sim_cmds <- data.frame(locus_number = locus_number,
                         command = args,
                         stringsAsFactors = FALSE,
                         row.names = NULL)

  reduce_sim_commands(sim_cmds)
}


reduce_sim_commands <- function(sim_commands) {
  if (nrow(sim_commands) == 1) return(sim_commands)
  grouped_commands <- unique(sim_commands[, 2])
  if (length(grouped_commands) == nrow(sim_commands)) return(sim_commands)
  grouped_locus_number <- vapply(grouped_commands, function(cmd) {
    sum(sim_commands[sim_commands[, 2] == cmd, 1])
  }, numeric(1)) #nolint
  data.frame(locus_number = grouped_locus_number,
             command = grouped_commands,
             stringsAsFactors = FALSE,
             row.names = NULL)
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
