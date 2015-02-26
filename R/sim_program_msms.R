# --------------------------------------------------------------
# sim_prog_msms.R
# Calling msms from a demographic model.
#
# Authors:  Paul R. Staab
# Date:     2014-06-30
# Licence:  GPLv3 or later
# --------------------------------------------------------------

msms.features <- c("selection", "selection_AA", "selection_Aa")
possible.features  <- c(get_sim_prog('ms')$possible_features, msms.features)
possible.sum.stats <- get_sim_prog('ms')$possible_sum_stats

callMsms <- function(ms.args, msms.args, subgroup) {
  checkForMsms()

  out_file <- tempfile('msms')
  seed <- sampleSeed(1)

  # Create the command
  cmd <- paste("java -jar", get_msms_path(), as.character(msms.args),
               "-ms", as.character(ms.args), "-seed", seed, ">", out_file)

  # Execute the command
  capture.output(system(cmd))

  if(!file.exists(out_file)) stop("msms simulation failed!")
  if(file.info(out_file)$size == 0) stop("msms output is empty!")

  out_file
}

checkForMsms <- function(throw.error = TRUE, silent = FALSE) {
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
generateMsmsOptionsCommand <- function(dm) {
  cmd <- c('c(')

  for (i in 1:dim(dm$features)[1] ) {
    type <- as.character(dm$features[i,"type"])
    feat <- unlist(dm$features[i, ])

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


generateMsmsOptions <- function(dm, parameters, locus) {
  msms.tmp <- createParameterEnv(dm, parameters, locus = locus)

  cmd <- read_cache(dm, 'msms_cmd')
  if (is.null(cmd)) {
    cmd <- generateMsmsOptionsCommand(dm)
    cache(dm, 'msms_cmd', cmd)
  }

  eval(parse(text=cmd), envir=msms.tmp)
}

msmsSimFunc <- function(dm, parameters) {
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
                        paste(generateMsOptions(dm, parameters, locus),
                              collapse=" "))
    msms.options <- paste(generateMsmsOptions(dm, parameters, locus),
                          collapse= " ")
    #print(c(ms.options, msms.options))
    callMsms(ms.options, msms.options)
  })

  # Parse the output and calculate summary statistics
  seg_sites <- parseMsOutput(files,
                             get_sample_size(dm, for_sim = TRUE),
                             get_locus_number(dm))

  sum_stats <- calc_sumstats(seg_sites, files, dm, parameters)

  # Clean Up
  unlink(files)
  sum_stats
}


msms_get_command <- function(model) {
  ms_cmd <- paste(generateMsOptions(model, get_parameter_table(model)$name),
                  collapse = ' ')
  msms_cmd <- paste(generateMsmsOptions(model,
                                        get_parameter_table(model)$name, 1),
                    collapse = ' ')

  paste("msms", msms_cmd,
        "-ms", sum(get_sample_size(model)), get_locus_number(model), ms_cmd)
}


#' @include sim_program.R
createSimProgram("msms", possible.features, possible.sum.stats,
                 msmsSimFunc, msms_get_command, priority=40)
