create_par_env <- function(model, parameters, ...) {
  par_env <- new.env()

  par_names <- get_parameter_table(model)$name
  for (i in seq(along = par_names)){
    par_env[[par_names[i]]] <- parameters[i]
  }

  additional_pars <- list(...)
  for (i in seq(along = additional_pars)) {
    par_env[[names(additional_pars)[i]]] <- additional_pars[[i]]
  }

  par_env[["locus_number"]] <- get_locus_number(model)

  par_env
}


escape_par_expr <- function(cmd) {
  if (is.null(names(cmd))) return(cmd)
  for (i in seq(along = cmd)) {
    if (names(cmd)[i] == '') next
    cmd[i] <- paste0('\"', cmd[i], '\"')
  }
  cmd
}
