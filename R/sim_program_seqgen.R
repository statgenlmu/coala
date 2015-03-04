# --------------------------------------------------------------
# simprog_seqgen.R
# Adaptor to calling ms from a demographic model.
#
# Authors:  Paul R. Staab
# Date:     2013-11-21
# Licence:  GPLv3 or later
# --------------------------------------------------------------

# list ms's features + FS related features
sg_features <- unique(c(#get_simprog('ms')$possible_features,
                        #get_simprog('msms')$possible_features,
                        'mutation_model', 'tstv_ratio',
                        'base_freq_A', 'base_freq_C', 'base_freq_G',
                        'base_freq_T',
                        'gtr_rate_1', 'gtr_rate_2', 'gtr_rate_3',
                        'gtr_rate_4','gtr_rate_5','gtr_rate_6',
                        'gamma_categories', 'gamma_rate',
                        'locus_trios', 'outgroup',
                        'mutation_outer'))

sg_sum_stats <- c('jsfs', 'file', 'seg.sites')
sg_mutation_models <- c('HKY', 'F84', 'GTR')

sg_find_exe <- function(throw.error = TRUE, silent = FALSE) {
  if ((!is.null(get_seqgen_path())) && file.exists(get_seqgen_path())) {
    return(TRUE)
  }

  # Works on Linux only maybe
  run.path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
  executables <- c(file.path(run.path, "seq-gen"),
                   file.path(run.path, "seqgen"))
  for (exe in executables) {
    if (file.exists(exe)) {
      if (!silent) message(paste("Using", exe, "as seqgen executable\n"))
      set_seqgen_path(exe)
      return(TRUE)
    }
  }

  if (throw.error) {
    stop("No seqgen executable found. Please provide one using
          set_seqgen_exe()")
  }
  return(FALSE)
}

generate_tree_model <- function(dm, locus, locus_number=1) {
  tree_model <- read_cache(dm, paste0('tree_model_', locus))

  if (is.null(tree_model)) {
    locus_length <- get_locus_length_matrix(dm)[locus,]

    if (any(msms_features %in% get_feature_table(dm)$type)) {
      tree.prog <- get_simprog('msms')
    } else {
      tree.prog <- get_simprog('ms')
    }

    tree_model <- dm

    # Features
    tree_model$features <-
      dm$features[get_feature_table(dm)$type %in% tree.prog$possible_features]

    # Summary Stastics
    tree_model$sum_stats <- create_sumstat_container()
    tree_model <- tree_model + sumstat_sg_trees(locus_length)

    # Loci
    tree_model$loci <- tree_model$loci[FALSE, ]
    if (locus_number == 1) {
      tree_model <- tree_model + locus_single(sum(locus_length))
    } else {
      tree_model <- tree_model + locus_averaged(locus_number,
                                                sum(locus_length))
    }

    cache(dm, paste0('tree_model_', locus), tree_model)
  }

  tree_model
}


#' Set the path to the executable for seqgen
#'
#' @param seqgen.exe Path to seqgen's executable.
#' @export
set_seqgen_exe <- function(seqgen_exe) {
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
sg_call <- function(opts, ms_files) {
  sg_find_exe()
  stopifnot(!missing(opts))
  stopifnot(length(opts) == length(ms_files))

  sapply(seq(along = opts), function(i) {
    if(file.info(ms_files[i])$size == 0 ) stop("No trees in file ", ms_files[i])

    seqgen_file <- tempfile('seqgen')
    cmd <- paste(opts[i], "<", ms_files[i], ">", seqgen_file)

    # Do the acctual simulation
    if (system(cmd, intern = FALSE) != 0) stop("seq-gen simulation failed")

    if( !file.exists(seqgen_file) ) stop("seq-gen simulation failed!")
    if( file.info(seqgen_file)$size == 0 ) stop("seq-gen output is empty!")

    seqgen_file
  })
}

sg_generate_opts <- function(dm, parameters, locus,
                                  seeds, eval_pars = TRUE) {

  cmd <- read_cache(dm, 'seqgen_cmd')
  if (is.null(cmd)) {
    cmd <- sg_generate_opt_cmd(dm)
    cache(dm, 'seqgen_cmd', cmd)
  }

  locus_lengths <- get_locus_length_matrix(dm)[locus,]
  if (locus_lengths[1] == 0 & locus_lengths[5] == 0) {
    locus_lengths <- locus_lengths[3]
  } else {
    locus_lengths <- locus_lengths[c(1,3,5)]
  }

  # Fill the parameters in the template
  sapply(seq(along = locus_lengths), function(i) {
    if (!eval_pars) cmd[[i]] <- escape_par_expr(cmd[[i]])
    par_envir <- create_par_env(dm, parameters, locus = locus,
                                    locus_length = locus_lengths[i],
                                    seed = seeds[i])
    paste(eval(parse(text=cmd[[i]]), envir=par_envir), collapse=" ")
  })
}


sg_generate_opt_cmd <- function(dm) {
  stopifnot(is.model(dm))

  if (length(get_outgroup_size(dm)) == 0) {
    stop("Finite Sites models need an outgroup.")
  }

  if (has_trios(dm)) is_outer <- c(TRUE, FALSE, TRUE)
  else is_outer <- FALSE

  lapply(is_outer, function(outer) {
    opts <- c('c(', paste('"', get_seqgen_path(), '"', sep=""), ",")
    base.freqs <- list()
    gtr.rates <- list()

    for (i in 1:dim(get_feature_table(dm))[1] ) {
      type <- as.character(get_feature_table(dm)[i,"type"])
      feat <- unlist(get_feature_table(dm)[i, ])

      if (type == "mutation_model") {
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

    opts <- c(opts, '"-l"', ',', 'locus_length', ',')
    opts <- c(opts, '"-s"', ',',
              s=paste(get_mutation_par(dm, outer), ' / locus_length'), ',')
    opts <- c(opts, '"-p"', ',', 'locus_length + 1', ',')
    opts <- c(opts, '"-z"', ',', 'seed', ',')
    opts <- c(opts, '"-q"', ')')
    opts
  })
}


sg_get_command <- function(dm) {
  tree_model <- generate_tree_model(dm, 1)
  tree_cmd <-
    get_simprog(determine_simprog(tree_model))$print_cmd_func(tree_model)

  sg_cmd <- paste(sg_generate_opts(dm, get_parameter_table(dm)$name,
                                        locus = 1, seeds = 'seed',
                                        eval_pars=FALSE),
                  collapse = ' ')
  sg_cmd <- paste('seq-gen', sg_cmd)

  c(tree=tree_cmd, seqgen=sg_cmd)
}


sg_simulate <- function(dm, parameters) {
  sg_find_exe()
  if (length(get_outgroup_size(dm)) == 0) {
    stop("Finite site models need an outgroup")
  }

  # Run all simulation in with one seqgen call if they loci are identical,
  # or call ms for each locus if there is variation between the loci.

  if (has_inter_locus_var(dm)) {
    sim_reps <- 1:get_locus_number(dm)
    sim_loci <- 1
  } else {
    sim_reps <- 1
    sim_loci <- get_locus_number(dm)
  }

  seqgen.files <- lapply(sim_reps, function(locus) {
    # Generate options for seqgen
    tree.model <- generate_tree_model(dm, locus, sim_loci)
    stopifnot(!is.null(tree.model))

    # Simulate the trees
    sum_stats_ms <- simulate(tree.model, pars=parameters)

    # Call seq-gen to distribute mutations
    seqgen.options <-
      sg_generate_opts(dm, parameters, locus,
                            sample_seed(length(sum_stats_ms$sg_trees)))

    seqgen.file <- sg_call(seqgen.options, sum_stats_ms$sg_trees)

    # Delete tree files
    unlink(c(sum_stats_ms[['file']], sum_stats_ms$sg_trees))
    seqgen.file
  })
  stopifnot(length(seqgen.files) == length(sim_reps))

  # Generate the summary statistics
  seg_sites <- parse_sg_output(seqgen.files,
                               sum(get_sample_size(dm, for_sim = TRUE)),
                               get_locus_length_matrix(dm),
                               get_locus_number(dm),
                               outgroup_size = get_outgroup_size(dm, TRUE))

  sum_stats <- calc_sumstats(seg_sites, seqgen.files, dm, parameters)

  # Clean Up
  unlink(unlist(seqgen.files))
  sum_stats
}


# @include sim_program.R
#create_simprog("seq-gen", sg_features, sg_sum_stats,
#                 sg_simulate, sg_get_command,
#                 priority=10)
