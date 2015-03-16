call_msms <- function(ms_args, msms_args, subgroup) {
  msms_find_jar()

  out_file <- tempfile('msms')
  seed <- sample_seed(1)

  # Create the command
  cmd <- paste("java -jar", get_msms_path(), as.character(msms_args),
               "-ms", as.character(ms_args), "-seed", seed, ">", out_file)

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



# This function generates an string that contains an R command for generating
# an msms call to the current model_
msms_generate_opts_cmd <- function(model) {
  cmd <- c('c(')

  for (i in 1:dim(get_feature_table(model))[1] ) {
    type <- as.character(get_feature_table(model)[i,"type"])
    feat <- unlist(get_feature_table(model)[i, ])

    if (type == "selection") {
      cmd <- c(cmd, '"-SI"', ',', feat['time.point'], ',',
               length(get_sample_size(model, for_sim = TRUE)), ',')
      start_freq <- rep(0, length(get_sample_size(model, for_sim = TRUE)))
      start_freq[ as.integer(feat['pop.source']) ] <- 0.0005
      cmd <- c(cmd, paste0('"', paste(start_freq, collapse=' '), '"'), ',')

      cmd <- c(cmd, '"-N 10000"', ',')

      s_AA <- search_feature(model, 'selection_AA',
                             pop.source = feat['pop.source'],
                             time.point = feat['time.point'])$parameter
      stopifnot(length(s_AA) == 1)
      cmd <- c(cmd, '"-SAA"', ',', s_AA, ',')

      s_Aa <- search_feature(model, 'selection_Aa',
                            pop.source = feat['pop.source'],
                            time.point = feat['time.point'])$parameter
      stopifnot(length(s_Aa) == 1)
      cmd <- c(cmd, '"-SAa"', ',', s_Aa, ',')
      cmd <- c(cmd, '"-Sp 0.5"', ',', '"-SForceKeep"', ',')
      cmd <- c(cmd, '"-threads 1"', ',')
    }
  }

  cmd <- c(cmd, '" ")')
  cmd
}


msms_generate_opts <- function(model, parameters, locus) {
  msms_tmp <- create_par_env(model, parameters,
                             locus = locus)

  cmd <- read_cache(model, 'msms_cmd')
  if (is.null(cmd)) {
    cmd <- msms_generate_opts_cmd(model)
    cache(model, 'msms_cmd', cmd)
  }

  eval(parse(text=cmd), envir=msms_tmp)
}


#' @include simulator_class.R
#' @include simulator_ms.R
Simulator_msms <- R6Class("Simulator_msms", inherit = Simulator,
  private = list(
    name = 'msms',
    features = c("selection", "selection_AA", "selection_Aa",
                 get_simulator("ms")$get_features()),
    sumstats = get_simulator("ms")$get_sumstats(),
    priority = 40
  ),
  public = list(
    get_cmd = function(model) {
      ms_cmd <- paste(ms_generate_opts(model,
                                       get_parameter_table(model)$name,
                                       "locus",
                                       "locus_length"),
                      collapse = ' ')
      msms_cmd <- paste(msms_generate_opts(model,
                                           get_parameter_table(model)$name,
                                           "locus"),
                        collapse = ' ')

      paste("msms", msms_cmd,
            "-ms", sum(get_sample_size(model)), get_locus_number(model), ms_cmd)
    },
    simulate = function(model, parameters) {
      # Get the length and number of loci
      llm <- get_locus_length_matrix(model, has_inter_locus_var(model))

      # Run the simulation(s)
      files <- lapply(1:nrow(llm), function(i) {
        ms_options <- paste(sum(get_sample_size(model, for_sim = TRUE)),
                            llm[i, 'number'],
                            paste(ms_generate_opts(model, parameters, i,
                                                   sum(llm[i, 1:5])),
                                  collapse=" "))
        msms_options <- paste(msms_generate_opts(model, parameters, i),
                              collapse= " ")
        #print(c(ms_options, msms_options))
        call_msms(ms_options, msms_options)
      })

      # Parse the output and calculate summary statistics
      if (requires_segsites(model)) {
        seg_sites <- parse_ms_output(files,
                                     get_sample_size(model, for_sim = TRUE),
                                     get_locus_number(model))
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
