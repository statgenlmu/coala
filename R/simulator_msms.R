call_msms <- function(ms.args, msms.args, subgroup) {
  msms_find_jar()

  out_file <- tempfile('msms')
  seed <- sample_seed(1)

  # Create the command
  cmd <- paste("java -jar", get_msms_path(), as.character(msms.args),
               "-ms", as.character(ms.args), "-seed", seed, ">", out_file)

  # Execute the command
  capture.output(system(cmd))

  if(!file.exists(out_file)) stop("msms simulation failed!")
  if(file.info(out_file)$size == 0) stop("msms output is empty!")

  out_file
}

msms_find_jar <- function(throw.error = TRUE, silent = FALSE) {
  if ((!is.null(get_msms_path())) && file.exists(get_msms_path())) return(TRUE)

  # Works on Linux only maybe
  run.path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
  executables <- file.path(c(run.path, getwd()), "msms.jar")
  for (exe in executables) {
    if (file.exists(exe)) {
      if (!silent) message(paste("Using", exe, "as msms implementation\n"))
      set_msms_path(exe)
      return(TRUE)
    }
  }

  if (throw.error) stop("No msms executable found.")
  FALSE
}



# This function generates an string that contains an R command for generating
# an msms call to the current model.
msms_generate_opts_cmd <- function(dm) {
  cmd <- c('c(')

  for (i in 1:dim(get_feature_table(dm))[1] ) {
    type <- as.character(get_feature_table(dm)[i,"type"])
    feat <- unlist(get_feature_table(dm)[i, ])

    if (type == "selection") {
      cmd <- c(cmd, '"-SI"', ',', feat['time.point'], ',',
               length(get_sample_size(dm, for_sim = TRUE)), ',')
      start.freq <- rep(0, length(get_sample_size(dm, for_sim = TRUE)))
      start.freq[ as.integer(feat['pop.source']) ] <- 0.0005
      cmd <- c(cmd, paste0('"', paste(start.freq, collapse=' '), '"'), ',')

      cmd <- c(cmd, '"-N 10000"', ',')

      s_AA <- search_feature(dm, 'selection_AA',
                            pop.source = feat['pop.source'],
                            time.point = feat['time.point'])$parameter
      stopifnot(length(s_AA) == 1)
      cmd <- c(cmd, '"-SAA"', ',', s_AA, ',')

      s_Aa <- search_feature(dm, 'selection_Aa',
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


msms_generate_opts <- function(dm, parameters, locus) {
  msms.tmp <- create_par_env(dm, parameters, locus = locus)

  cmd <- read_cache(dm, 'msms_cmd')
  if (is.null(cmd)) {
    cmd <- msms_generate_opts_cmd(dm)
    cache(dm, 'msms_cmd', cmd)
  }

  eval(parse(text=cmd), envir=msms.tmp)
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
      ms_cmd <- paste(ms_generate_opts(model, get_parameter_table(model)$name),
                      collapse = ' ')
      msms_cmd <- paste(msms_generate_opts(model,
                                           get_parameter_table(model)$name, 1),
                        collapse = ' ')

      paste("msms", msms_cmd,
            "-ms", sum(get_sample_size(model)), get_locus_number(model), ms_cmd)
    },
    simulate = function(dm, parameters) {
      # Run all simulation in with one ms call if they loci are identical,
      # or call ms for each locus if there is variation between the loci.
      if (has_inter_locus_var(dm)) {
        sim_reps <- 1:get_locus_number(dm)
        sim_loci <- 1
      } else {
        sim_reps <- 1
        sim_loci <- get_locus_number(dm)
      }

      # Run the simulation(s)
      files <- lapply(sim_reps, function(locus) {
        ms.options <- paste(sum(get_sample_size(dm, for_sim = TRUE)), sim_loci,
                            paste(ms_generate_opts(dm, parameters, locus),
                                  collapse=" "))
        msms.options <- paste(msms_generate_opts(dm, parameters, locus),
                              collapse= " ")
        #print(c(ms.options, msms.options))
        call_msms(ms.options, msms.options)
      })

      # Parse the output and calculate summary statistics
      seg_sites <- parse_ms_output(files,
                                   get_sample_size(dm, for_sim = TRUE),
                                   get_locus_number(dm))

      sum_stats <- calc_sumstats(seg_sites, files, dm, parameters)

      # Clean Up
      unlink(files)
      sum_stats
    }
  )
)

register_simulator(Simulator_msms)

