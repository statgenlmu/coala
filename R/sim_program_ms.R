# --------------------------------------------------------------
# Translates an demographic model to an ms command and 
# executes the simulation.
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Licence:  GPLv3 or later
# --------------------------------------------------------------

possible.features  <- c("sample", "mutation", "migration", "migration_sym",
                        "pop_merge",
                        "recombination", "size_change", "growth",
                        "inter_locus_variation")
possible.sum.stats <- c("jsfs", "trees", "seg.sites", "file")


# This function generates an string that contains an R command for generating
# an ms call to the current model.
generateMsOptionsCommand <- function(dm) {
  nSample <- dm.getSampleSize(dm)
  cmd <- c('c(')
  cmd <- c(cmd,'"-I"', ",", length(nSample), ',', 
           paste(nSample, collapse=","), ',')

  for (i in 1:dim(dm@features)[1] ) {
    type <- as.character(dm@features[i,"type"])
    feat <- unlist(dm@features[i, ])

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
      cmd <- c(cmd, '"-r"', ',', feat['parameter'], ',', dm.getLociLength(dm), ',')

    else if (type == "size.change"){
      cmd <- c(cmd, '"-en"', ',', feat['time.point'], ',',
               feat["pop.source"], ',', feat['parameter'], ',')
    }

    else if (type == "growth"){
      cmd <- c(cmd, '"-eg"', ',' , feat["time.point"], ',',
               feat["pop.source"], ',', feat["parameter"], ',')
      }

    else if (type %in% c("sample", "loci.number", "loci.length", 
                         "pos.selection", "bal.selection",
                         "inter_locus_variation")) {}
    else stop("Unknown feature:", type)
  }

  if ('trees' %in% dm.getSummaryStatistics(dm)) cmd <- c(cmd, '"-T",')
  cmd <- c(cmd, '" ")')
}

generateMsOptions <- function(dm, parameters, subgroup) {
  ms.tmp <- new.env()

  par.names <- dm.getParameters(dm)
  for (i in seq(along = par.names)){
    ms.tmp[[ par.names[i] ]] <- parameters[i]
  }

  fixed.pars <- dm@parameters[dm@parameters$fixed, ]
  if (nrow(fixed.pars) > 0) {
    for (i in 1:nrow(fixed.pars)){
      ms.tmp[[ fixed.pars$name[i] ]] <- fixed.pars$lower.range[i]
    }
  }

  if ( !is.null( dm@options[['ms.cmd']] ) )
    cmd <- dm@options[['ms.cmd']]
  else
    cmd <- generateMsOptionsCommand(dm)
  cmd <- eval(parse(text=cmd), envir=ms.tmp)

  return(cmd)
}

printMsCommand <- function(dm) {
  cmd <- generateMsOptionsCommand(dm)

  cmd <- cmd[cmd != ","]
  cmd <- cmd[-c(1, length(cmd))]

  cmd <- paste(cmd, collapse=" ")

  cmd <- gsub(",", " ", cmd)
  cmd <- gsub('\"', "", cmd)
  cmd <- gsub('"', " ", cmd)

  cmd <- paste("ms", sum(dm.getSampleSize(dm)), dm.getLociNumber(dm), cmd)
  .print(cmd)
}

msSingleSimFunc <- function(dm, parameters) {
  checkType(dm, "dm")
  checkType(parameters, "num")
  if (length(parameters) != dm.getNPar(dm)) stop("Wrong number of parameters!")

  # Run all simulation in with one ms call if they loci are identical,
  # or call ms for each locus if there is variation between the loci.
  if (hasInterLocusVariation(dm)) {
    sim_reps <- 1:dm.getLociNumber(dm)
    sim_loci <- 1
  } else {
    sim_reps <- 1
    sim_loci <- dm.getLociNumber(dm)
  }
  
  # Do the actuall simulation
  ms.files <- lapply(sim_reps, function(locus) {
    ms.options <- generateMsOptions(dm, parameters, locus)
    ms.file <- getTempFile('ms')
    ms(sum(dm.getSampleSize(dm)), sim_loci, 
       unlist(strsplit(ms.options, " ")), ms.file)
    ms.file
  })
  
  # Parse & return the simulation output
  generateSumStats(ms.files, 'ms', parameters, dm)
}

finalizeMs <- function(dm) {
  dm@options[['ms.cmd']] <- generateMsOptionsCommand(dm)
  return(dm)
}

#' @include dm_sim_program.R
createSimProgram("ms", possible.features, possible.sum.stats, 
                 msSingleSimFunc, finalizeMs, printMsCommand, 100)