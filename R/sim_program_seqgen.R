# --------------------------------------------------------------
# sim_prog_seqgen.R
# Adaptor to calling ms from a demographic model.
#
# Authors:  Paul R. Staab
# Date:     2013-11-21
# Licence:  GPLv3 or later
# --------------------------------------------------------------

# list ms's features + FS related features
sg.features <- unique(c(getSimProgram('ms')$possible_features,
                        getSimProgram('msms')$possible_features,
                        'mutation_model', 'tstv_ratio',
                        'base_freq_A', 'base_freq_C', 'base_freq_G',
                        'base_freq_T',
                        'gtr_rate_1', 'gtr_rate_2', 'gtr_rate_3',
                        'gtr_rate_4','gtr_rate_5','gtr_rate_6',
                        'gamma_categories', 'gamma_rate',
                        'locus_trios', 'outgroup',
                        'mutation_outer'))

sg.sum.stats <- c('jsfs', 'file', 'seg.sites')
sg.mutation.models <- c('HKY', 'F84', 'GTR')

checkForSeqgen <- function(throw.error = TRUE, silent = FALSE) {
  if ((!is.null(get_seqgen_path())) && file.exists(get_seqgen_path())) {
    return(TRUE)
  }

  # Works on Linux only maybe
  run.path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
  executables <- c(paste(run.path, "/seq-gen", sep=""),
                   paste(run.path, "/seqgen", sep=""))
  for (exe in executables) {
    if (file.exists(exe)) {
      if (!silent) message(paste("Using", exe, "as seqgen executable\n"))
      set_seqgen_path(exe)
      return(TRUE)
    }
  }

  if (throw.error) {
    stop("No seqgen executable found. Please provide one using
          setSeqgenExecutable()")
  }
  return(FALSE)
}

generateTreeModel <- function(dm, locus_length) {
  stopifnot(all(get_groups(dm) == 1))
  if (any(msms.features %in% dm$features$type)) {
    tree.prog <- getSimProgram('msms')
  } else {
    tree.prog <- getSimProgram('ms')
  }

  dm$features <- dm$features[dm$features$type %in% tree.prog$possible_features, ]
  dm <- resetSumStats(dm)
  dm <- dm.addSummaryStatistic(dm, "trees")
  dm <- dm.addSummaryStatistic(dm, "file")

  dm$loci <- dm$loci[FALSE, ]
  dm <- addLocus(dm, number = 1, length_m = sum(locus_length), group = 0)
  dm.finalize(dm)
}


#' Set the path to the executable for seqgen
#'
#' @param seqgen.exe Path to seqgen's executable.
#' @export
setSeqgenExecutable <- function(seqgen_exe) {
  if (file.exists(seqgen_exe)) {
    set_seqgen_path(seqgen_exe)
  } else {
    stop("File", seqgen_exe, "does not exist")
  }
}

# Function to perform simulation using seqgen
#
# @param opts The options to pass to ms. Must either be a character or character
# vector.
callSeqgen <- function(opts, ms_files) {
  stopifnot(!missing(opts))
  stopifnot(length(opts) == length(ms_files))

  sapply(seq(along = opts), function(i) {
    if(!file.exists(ms_files[i])) stop("ms file not found")
    if(file.info(ms_files[i])$size == 0 ) stop("ms output is empty")

    seqgen_file <- tempfile('csr_seqgen')
    cmd <- paste(opts[i], "<", ms_files[i], ">", seqgen_file)

    # Do the acctual simulation
    if (system(cmd, intern = FALSE) != 0) stop("seq-gen simulation failed")

    if( !file.exists(seqgen_file) ) stop("seq-gen simulation failed!")
    if( file.info(seqgen_file)$size == 0 ) stop("seq-gen output is empty!")

    seqgen_file
  })
}

generateSeqgenOptions <- function(dm, parameters, locus,
                                  locus_lengths, seeds) {
  # Generate the command template to execute or use the buffered one
  if ( !is.null( dm$options[['seqgen.cmd']] ) ) {
    cmd <- dm$options[['seqgen.cmd']]
  } else {
    cmd <- generateSeqgenOptionsCmd(dm)
  }

  if (locus_lengths[1] == 0 & locus_lengths[5] == 0) {
    locus_lengths <- locus_lengths[3]
  } else {
    locus_lengths <- locus_lengths[c(1,3,5)]
  }

  # Fill the parameters in the template
  sapply(seq(along = locus_lengths), function(i) {
    par_envir <- createParameterEnv(dm, parameters, locus = locus,
                                    locus_length = locus_lengths[i],
                                    seed = seeds[i])
    paste(eval(parse(text=cmd[[i]]), envir=par_envir), collapse=" ")
  })
}


generateSeqgenOptionsCmd <- function(dm) {
  stopifnot(is.model(dm))
  base.freqs <- F
  gtr.rates <- F
  includes.model <- F

  if (!is.numeric(get_outgroup_size(dm))) {
    stop("Finite Sites models need an outgroup.")
  }

  if (!dm.hasTrios(dm)) is_outer <- FALSE
  else is_outer <- c(TRUE, FALSE, TRUE)

  lapply(is_outer, function(outer) {
    opts <- c('c(', paste('"', get_seqgen_path(), '"', sep=""), ",")
    base.freqs <- list()
    gtr.rates <- list()

    for (i in 1:dim(dm$features)[1] ) {
      type <- as.character(dm$features[i,"type"])
      feat <- unlist(dm$features[i, ])

      if (type == "mutation_model") {
        includes.model <- T
        opts <- c(opts, paste('"-m', feat['parameter'], '"', sep=""), ",")
      }

      else if ( type %in% c('base_freq_A', 'base_freq_C',
                            'base_freq_G', 'base_freq_T') )
        base.freqs[[type]] <- feat['parameter']

      else if ( type %in% c('gtr_rate_1', 'gtr_rate_2', 'gtr_rate_3',
                            'gtr_rate_4', 'gtr_rate_5', 'gtr_rate_6') )
        gtr.rates[[type]] <- feat['parameter']

      else if (type == "tstv_ratio")
        opts <- c(opts, '"-t"', ',', feat['parameter'], ',')

      else if (type == "gamma_rate")
        opts <- c(opts, '"-a"', ',', feat['parameter'], ',')

      else if (type == "gamma_categories")
        opts <- c(opts, '"-g"', ',', feat['parameter'], ',')
    }

    if (length(base.freqs) == 4) {
      opts <- c(opts, '"-f"', ',', base.freqs[['base_freq_A']],
                ',', base.freqs[['base_freq_C']],
                ',', base.freqs[['base_freq_G']],
                ',', base.freqs[['base_freq_T']], ',')
    }

    if (length(gtr.rates) == 6) {
      opts <- c(opts, '"-r"', ',', gtr.rates[['gtr_rate_1']],
                ',', gtr.rates[['gtr_rate_2']],
                ',', gtr.rates[['gtr_rate_3']],
                ',', gtr.rates[['gtr_rate_4']],
                ',', gtr.rates[['gtr_rate_5']],
                ',', gtr.rates[['gtr_rate_6']], ',')
    }

    if (!includes.model) {
      stop("You must specify a finite sites mutation model for this demographic model")
    }

    opts <- c(opts, '"-l"', ',', 'locus_length', ',')
    opts <- c(opts, '"-s"', ',', paste(getThetaName(dm, outer), ' / locus_length'), ',')
    opts <- c(opts, '"-p"', ',', 'locus_length + 1', ',')
    opts <- c(opts, '"-z"', ',', 'seed', ',')
    opts <- c(opts, '"-q"', ')')
    opts
  })
}

printSeqgenCommand <- function(dm) {
  tree.model <- generateTreeModel(dm, get_locus_length_matrix(dm)[1,3])
  getSimProgram(tree.model$currentSimProg)$print_cmd_func(tree.model)

  cmds <- generateSeqgenOptionsCmd(dm)
  for (cmd in cmds) {
    cmd <- cmd[cmd != ","]
    cmd <- cmd[-c(1, length(cmd))]

    cmd <- paste(cmd, collapse=" ")

    cmd <- gsub(",", " ", cmd)
    cmd <- gsub('\"', "", cmd)
    cmd <- gsub('"', " ", cmd)

    cat(cmd, '\n')
  }
}

seqgenSingleSimFunc <- function(dm, parameters) {
  checkForSeqgen()

  locus_length <- get_locus_length_matrix(dm)

  seqgen.files <- lapply(1:get_locus_number(dm), function(locus) {
    # Generate options for seqgen
    tree.model <- generateTreeModel(dm, locus_length[locus,])

    # Simulate the trees
    sum_stats_ms <- simulate(tree.model, parameters)
    tree_files <- parseTrees(sum_stats_ms[['file']][[1]],
                             locus_length[locus,],
                             tempfile)

    # Call seq-gen to distribute mutations
    seqgen.options <- generateSeqgenOptions(dm, parameters, locus,
                                            locus_length[locus,],
                                            sampleSeed(length(tree_files)))
    seqgen.file <- callSeqgen(seqgen.options, tree_files)

    # Delete tree files
    unlink(c(tree_files, sum_stats_ms[['file']]))
    seqgen.file
  })

  # Generate the summary statistics
  generateSumStats(seqgen.files, 'seqgen', parameters, dm)
}

finalizeSeqgen <- function(dm) {
  stopifnot(is.model(dm))
  checkForSeqgen()
  dm$options[['seqgen.cmd']] <- generateSeqgenOptionsCmd(dm)
  dm
}

#' @include sim_program.R
createSimProgram("seq-gen", sg.features, sg.sum.stats,
                 seqgenSingleSimFunc, finalizeSeqgen, printSeqgenCommand,
                 priority=10)
