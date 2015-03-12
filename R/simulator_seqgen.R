# --------------------------------------------------------------
# simprog_seqgen.R
# Adaptor to calling ms from a demographic model.
#
# Authors:  Paul R. Staab
# Date:     2013-11-21
# Licence:  GPLv3 or later
# --------------------------------------------------------------


sg_only_features <- c('mutation_model', 'tstv_ratio',
                      'base_freq_A', 'base_freq_C', 'base_freq_G',
                      'base_freq_T',
                      'gtr_rate_1', 'gtr_rate_2', 'gtr_rate_3',
                      'gtr_rate_4','gtr_rate_5','gtr_rate_6',
                      'gamma_categories', 'gamma_rate',
                      'locus_trios', 'outgroup',
                      'mutation_outer',
                      'sumstat_dna')

#' @include simulator_ms.R
#' @include simulator_msms.R
#' @include simulator_scrm.R
sg_features <- unique(c(get_simulator("ms")$get_features(),
                        get_simulator("msms")$get_features(),
                        get_simulator("scrm")$get_features(),
                        sg_only_features))

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
    stop("No seqgen executable found. Please provide one using ",
         "set_seqgen_exe()")
  }

  FALSE
}

generate_tree_model <- function(model) {
  tree_model <- read_cache(model, "tree_model")

  if (is.null(tree_model)) {
    tree_model <- model

    # Features
    tree_model$features <-
      model$features[!get_feature_table(model)$type %in% sg_only_features]

    # Summary Stastics
    tree_model$sum_stats <- create_sumstat_container()
    tree_model <- tree_model + sumstat_sg_trees()

    cache(model, "tree_model", tree_model)
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
  } else stop("File", seqgen_exe, "does not exist")
}

# Function to perform simulation using seqgen
#
# @param opts The options to pass to ms. Must either be a character or character
# vector.
sg_call <- function(opts, ms_files) {
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

sg_generate_opts <- function(model, parameters, locus,
                             locus_lengths, seeds, eval_pars = TRUE) {

  cmd <- read_cache(model, 'seqgen_cmd')
  if (is.null(cmd)) {
    cmd <- sg_generate_opt_cmd(model)
    cache(model, 'seqgen_cmd', cmd)
  }

  if (locus_lengths[1] == 0 & locus_lengths[5] == 0) {
    locus_lengths <- locus_lengths[3]
  } else {
    locus_lengths <- locus_lengths[c(1,3,5)]
  }

  # Fill the parameters in the template
  sapply(seq(along = locus_lengths), function(i) {
    if (!eval_pars) cmd[[i]] <- escape_par_expr(cmd[[i]])
    par_envir <- create_par_env(model, parameters, locus = locus,
                                    locus_length = locus_lengths[i],
                                    seed = seeds[i])
    paste(eval(parse(text=cmd[[i]]), envir=par_envir), collapse=" ")
  })
}


sg_generate_opt_cmd <- function(model) {
  stopifnot(is.model(model))

  if (has_trios(model)) is_outer <- c(TRUE, FALSE, TRUE)
  else is_outer <- FALSE

  lapply(is_outer, function(outer) {
    opts <- c('c(', paste('"', get_seqgen_path(), '"', sep=""), ",")
    base.freqs <- list()
    gtr.rates <- list()

    for (i in 1:dim(get_feature_table(model))[1] ) {
      type <- as.character(get_feature_table(model)[i,"type"])
      feat <- unlist(get_feature_table(model)[i, ])

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
              s=paste(get_mutation_par(model, outer), ' / locus_length'), ',')
    opts <- c(opts, '"-p"', ',', 'locus_length + 1', ',')
    opts <- c(opts, '"-z"', ',', 'seed', ',')
    opts <- c(opts, '"-q"', ')')
    opts
  })
}


sg_get_command <- function(model) {
  tree_model <- generate_tree_model(model)

  tree_cmd <-
    determine_simprog(tree_model)$get_cmd(tree_model)

  sg_cmd <- paste(sg_generate_opts(model, get_parameter_table(model)$name,
                                   get_locus_length_matrix(model)[1, 1:5],
                                   locus = 1, seeds = 'seed',
                                   eval_pars=FALSE),
                  collapse = ' ')
  sg_cmd <- paste('seq-gen', sg_cmd)

  c(tree=tree_cmd, seqgen=sg_cmd)
}


sg_simulate <- function(model, parameters) {
  sg_find_exe()

  # Simulate the ancestral trees
  tree_model <- generate_tree_model(model)
  trees <- simulate(tree_model, pars=parameters)$trees
  assert_that(!is.null(trees))

  # Get loci length and number
  llm <- get_locus_length_matrix(model, has_inter_locus_var(model))

  # Call seq-gen for each locus (trio)
  seqgen_files <- lapply(1:length(trees), function(locus) {
    seqgen_options <- sg_generate_opts(model, parameters, locus,
                                       llm[locus, 1:5],
                                       sample_seed(length(trees[[locus]])))
    sg_call(seqgen_options, trees[[locus]])
  })
  unlink(unlist(trees))

  # Generate the summary statistics
  if (any(get_summary_statistics(model) != "dna")) {
    seg_sites <- parse_sg_output(seqgen_files,
                                 sum(get_sample_size(model, for_sim = TRUE)),
                                 llm[, 1:5, drop=FALSE],
                                 get_locus_number(model),
                                 outgroup_size = get_outgroup_size(model, TRUE))
  } else {
    seg_sites <- NULL
  }

  sum_stats <- calc_sumstats(seg_sites, seqgen_files, model, parameters)

  # Clean Up
  unlink(unlist(seqgen_files))
  sum_stats
}


#' @importFrom R6 R6Class
#' @include simulator_class.R
Simulator_seqgen <- R6Class('Simulator_seqgen', inherit = Simulator,
  private = list(
    name = 'seqgen',
    features = sg_features,
    sumstats = sg_sum_stats,
    priority = 10
   ),
   public = list(
     simulate = sg_simulate,
     get_cmd = sg_get_command
   )
)


register_simulator(Simulator_seqgen)
