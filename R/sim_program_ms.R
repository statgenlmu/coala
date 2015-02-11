# --------------------------------------------------------------
# Translates an demographic model to an ms command and
# executes the simulation.
#
# Authors:  Lisha Mathew & Paul R. Staab
# Licence:  GPLv3 or later
# --------------------------------------------------------------

possible.features  <- c("sample", "mutation", "migration", "migration_sym",
                        "pop_merge", "recombination", "size_change", "growth",
                        "inter_locus_variation")
possible.sum.stats <- c("jsfs", "trees", "seg.sites", "file")


# This function generates an string that contains an R command for generating
# an ms call to the current model.
generateMsOptionsCommand <- function(dm) {
  nSample <- get_sample_size(dm)
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

    else if (type %in% c("sample", "loci.number", "loci.length",
                         "selection", "selection_AA", "selection_Aa",
                         "inter_locus_variation")) {}
    else stop("Unknown feature:", type)
  }

  if ('trees' %in% get_summary_statistics(dm)) cmd <- c(cmd, '"-T",')
  cmd <- c(cmd, '" ")')
}

generateMsOptions <- function(dm, parameters, subgroup) {
  ms.tmp <- new.env()

  par.names <- get_parameter_table(dm)$name
  for (i in seq(along = par.names)){
    ms.tmp[[ par.names[i] ]] <- parameters[i]
  }

  cmd <- read_cache(dm, 'ms_cmd')
  if (is.null(cmd)) {
    cmd <- generateMsOptionsCommand(dm)
    cache(dm, 'ms_cmd', cmd)
  }

  eval(parse(text=cmd), envir=ms.tmp)
}

printMsCommand <- function(dm) {
  cmd <- generateMsOptionsCommand(dm)

  cmd <- cmd[cmd != ","]
  cmd <- cmd[-c(1, length(cmd))]

  cmd <- paste(cmd, collapse=" ")

  cmd <- gsub(",", " ", cmd)
  cmd <- gsub('\"', "", cmd)
  cmd <- gsub('"', " ", cmd)

  cmd <- paste("ms", sum(get_sample_size(dm)), get_locus_number(dm), cmd)
  cat(cmd, "\n")
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
  ms.files <- lapply(sim_reps, function(locus) {
    ms.options <- generateMsOptions(dm, parameters, locus)
    ms.file <- tempfile('csr_ms')

    ms(sum(get_sample_size(dm)), sim_loci,
       unlist(strsplit(ms.options, " ")), ms.file)

    if(file.info(ms.file)$size == 0) stop("ms simulation output is empty")
    ms.file
  })

  # Parse & return the simulation output
  generateSumStats(ms.files, 'ms', parameters, dm)
}


#' @include sim_program.R
createSimProgram("ms", possible.features, possible.sum.stats,
                 msSingleSimFunc, printMsCommand, 100)
