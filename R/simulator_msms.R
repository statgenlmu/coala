call_msms <- function(msms_args) {
  msms_find_jar()

  out_file <- tempfile('msms')
  seed <- sample_seed(1)

  # Create the command
  cmd <- paste("java -jar", get_msms_path(), as.character(msms_args),
               "-seed", seed, ">", out_file)

  # Execute the command
  capture.output(system(cmd))

  if(!file.exists(out_file)) stop("msms simulation failed!")
  if(file.info(out_file)$size == 0) stop("msms output is empty!")

  out_file
}

msms_find_jar <- function(throw_error = TRUE, silent = FALSE) {
  if ((!is.null(get_msms_path())) && file.exists(get_msms_path())) return(TRUE)

  # Works on Linux only maybe
  run_path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
  executables <- file.path(c(run_path, getwd()), "msms.jar")
  for (exe in executables) {
    if (file.exists(exe)) {
      if (!silent) message(paste("Using", exe, "as msms implementation\n"))
      set_msms_path(exe)
      return(TRUE)
    }
  }

  if (throw_error) stop("No msms executable found_")
  FALSE
}


# Translating the model into simulation commands
conv_to_msms_arg <- function(feature, model) UseMethod("conv_to_msms_arg")
conv_to_msms_arg.default <- function(feature, model) {
  stop("Unknown feature when generating ms command")
}


# This function generates an string that contains an R command for generating
# an ms call to the current model.
msms_generate_opts_cmd <- function(model) {
  cmd <- paste(vapply(model$features, conv_to_msms_arg,
                      FUN.VALUE = character(1), model),
               collapse = "")
  paste0("c('",
         sum(get_sample_size(model, TRUE)),
         "', format(locus_number, scientific = FALSE), '", cmd, "')")
}


msms_generate_opts <- function(model, parameters, locus, eval_pars = TRUE) {
  msms_tmp <- create_par_env(model, parameters,
                             locus_length = get_locus_length(model,
                                                             group = locus),
                             locus_number = get_locus_number(model, locus),
                             locus = locus, for_cmd = !eval_pars)

  cmd <- read_cache(model, 'msms_cmd')
  if (is.null(cmd)) {
    cmd <- msms_generate_opts_cmd(model)
    cache(model, 'msms_cmd', cmd)
  }

  eval(parse(text = cmd), envir = msms_tmp)
}


#' @include simulator_class.R
#' @include simulator_ms.R
Simulator_msms <- R6Class("Simulator_msms", inherit = Simulator,
  private = list(
    name = 'msms',
    priority = 40
  ),
  public = list(
    get_cmd = function(model) {
      cmd <- paste(msms_generate_opts(model, NULL, 1, FALSE),
                   collapse = ' ')

      paste("msms", cmd)
    },
    simulate = function(model, parameters=numeric(0)) {
      # Run the simulation(s)
      files <- lapply(1:get_locus_group_number(model), function(i) {
        msms_options <- paste(msms_generate_opts(model, parameters, i),
                              collapse= " ")
        call_msms(msms_options)
      })

      # Parse the output and calculate summary statistics
      if (requires_segsites(model)) {
        seg_sites <- parse_ms_output(files,
                                     get_sample_size(model, for_sim = TRUE),
                                     get_locus_number(model))

        if (has_trios(model)) {
          seg_sites <- conv_for_trios(seg_sites, model)
        }
      } else {
        seg_sites <- NULL
      }

      sum_stats <- calc_sumstats(seg_sites, files, model, parameters)

      # Clean Up
      unlink(files)
      sum_stats
    }
  )
)

register_simulator(Simulator_msms)
