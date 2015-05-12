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


conv_to_msms_arg <- function(feature, model) UseMethod("conv_to_msms_arg")

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.default <- function(feature, model) {
  stop("Unknown feature when generating ms command")
}


msms_create_cmd_tempalte <- function(model) {
  cmd <- read_cache(model, 'msms_cmd')
  if (is.null(cmd)) {
    cmd <- paste(vapply(model$features, conv_to_msms_arg,
                        FUN.VALUE = character(1), model),
                 collapse = "")
    cmd <- paste0("c('", cmd, "')")
    cache(model, 'msms_cmd', cmd)
  }
  cmd
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
      template <- msms_create_cmd_tempalte(model)
      cmd <- fill_cmd_template(template, model, NULL, 1, eval_pars = FALSE)
      paste("msms",
            sum(get_sample_size(model, TRUE)),
            cmd[1, "locus_number"],
            cmd[1, "command"])
    },
    simulate = function(model, parameters=numeric(0)) {
      template <- msms_create_cmd_tempalte(model)
      sample_size <- sum(get_sample_size(model, for_sim = TRUE))

      # Run the simulation(s)
      files <- lapply(1:get_locus_group_number(model), function(i) {
        cmds <- fill_cmd_template(template, model, parameters, i)

        vapply(1:nrow(cmds), function(i) {
          msms_options <- paste(sample_size,
                                cmds[i, "locus_number"],
                                cmds[i, "command"])
          call_msms(msms_options)
        }, character(1))
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
      unlink(unlist(files))
      sum_stats
    }
  )
)

register_simulator(Simulator_msms)
