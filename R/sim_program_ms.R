# --------------------------------------------------------------
# Translates an demographic model to an ms command and
# executes the simulation.
#
# Authors:  Lisha Mathew & Paul R. Staab
# Licence:  GPLv3 or later
# --------------------------------------------------------------

possible.features  <- c("sample", "mutation", "migration", "migration_sym",
                        "pop_merge", "recombination", "size_change", "growth",
                        "inter_locus_variation", "trees",
                        "unphased", "ploidy", "samples_per_ind")
possible.sum.stats <- c("jsfs", "trees", "seg.sites", "file")


# This function generates an string that contains an R command for generating
# an ms call to the current model.
generateMsOptionsCommand <- function(dm) {
  nSample <- get_sample_size(dm, for_sim = TRUE)
  cmd <- c('c(')
  cmd <- c(cmd,'"-I"', ",", length(nSample), ',',
           paste(nSample, collapse=","), ',')

  for (i in 1:dim(dm$features)[1] ) {
    type <- as.character(dm$features[i,"type"])
    feat <- unlist(dm$features[i, ])

    if ( type == "mutation" ) {
      cmd <- c(cmd,'"-t"', ',', feat["parameter"], ',')
    }

    else if (type == "pop_merge") {
      cmd <- c(cmd, '"-ej"', ',', feat["time.point"], ',',
               feat["pop.source"], ',', feat["pop.sink"], ',')
    }

    else if (type == "migration")
      cmd <- c(cmd, '"-em"', ',', feat['time.point'], ',',
               feat['pop.sink'], ',', feat['pop.source']  , ',',
               feat['parameter'], ',')

    else if (type == "migration_sym")
      cmd <- c(cmd, '"-eM"', ',',
               feat['time.point'], ',',
               feat['parameter'], ',')

    else if (type == "recombination")
      cmd <- c(cmd, '"-r"', ',', feat['parameter'], ',', get_locus_length(dm), ',')

    else if (type == "size_change"){
      cmd <- c(cmd, '"-en"', ',', feat['time.point'], ',',
               feat["pop.source"], ',', feat['parameter'], ',')
    }

    else if (type == "growth"){
      cmd <- c(cmd, '"-eg"', ',' , feat["time.point"], ',',
               feat["pop.source"], ',', feat["parameter"], ',')
      }

    else if (type == 'trees') {
      cmd <- c(cmd, '"-T",')
    }

    else if (type %in% c("sample", "loci.number", "loci.length",
                         "selection", "selection_AA", "selection_Aa",
                         "inter_locus_variation", "unphased",
                         "ploidy", "samples_per_ind")) {}
    else stop("Unknown feature:", type)
  }


  cmd <- c(cmd, '" ")')
}

generateMsOptions <- function(dm, parameters, eval_pars = TRUE) {
  ms.tmp <- createParameterEnv(dm, parameters)

  cmd <- read_cache(dm, 'ms_cmd')
  if (is.null(cmd)) {
    cmd <- generateMsOptionsCommand(dm)
    cache(dm, 'ms_cmd', cmd)
  }

  if (!eval_pars) cmd <- escape_par_expr(cmd)
  eval(parse(text=cmd), envir=ms.tmp)
}


ms_get_command <- function(model) {
  cmd <- generateMsOptions(model, get_parameter_table(model)$name, FALSE)
  txt <- paste(cmd, collapse = ' ')
  paste("ms", sum(get_sample_size(model)), get_locus_number(model), txt)
}


#' @importFrom phyclust ms
msSingleSimFunc <- function(dm, parameters=numeric()) {
  stopifnot(length(parameters) == 0 | all(is.numeric(parameters)))

  # Run all simulation in with one ms call if they loci are identical,
  # or call ms for each locus if there is variation between the loci.
  if (hasInterLocusVariation(dm)) {
    sim_reps <- 1:get_locus_number(dm)
    sim_loci <- 1
  } else {
    sim_reps <- 1
    sim_loci <- get_locus_number(dm)
  }

  # Do the actuall simulation
  files <- lapply(sim_reps, function(locus) {
    ms.options <- generateMsOptions(dm, parameters, locus)
    file <- tempfile('csr_ms')

    ms(sum(get_sample_size(dm, for_sim = TRUE)), sim_loci,
       unlist(strsplit(ms.options, " ")), file)

    if(file.info(file)$size == 0) stop("ms simulation output is empty")
    file
  })

  # Parse the output and calculate summary statistics
  seg_sites <- parseMsOutput(files,
                             get_sample_size(dm, for_sim = TRUE),
                             get_locus_number(dm))

  if (is_unphased(dm)) seg_sites <- unphase_segsites(seg_sites,
                                                     get_ploidy(dm),
                                                     get_samples_per_ind(dm))

  sum_stats <- calc_sumstats(seg_sites, files, dm, parameters)

  # Clean Up
  unlink(files)
  sum_stats
}


#' @include sim_program.R
createSimProgram("ms", possible.features, possible.sum.stats,
                 msSingleSimFunc, ms_get_command, 100)
