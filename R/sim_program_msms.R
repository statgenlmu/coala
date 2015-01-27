# --------------------------------------------------------------
# sim_prog_msms.R
# Calling msms from a demographic model.
# 
# Authors:  Paul R. Staab
# Date:     2014-06-30
# Licence:  GPLv3 or later
# --------------------------------------------------------------

msms.features <- c("pos.selection", "bal.selection")
possible.features  <- c(getSimProgram('ms')$possible_features, msms.features)
possible.sum.stats <- getSimProgram('ms')$possible_sum_stats

callMsms <- function(jar.path, ms.args, msms.args, subgroup) {
  out.file = getTempFile("msms")
  seed <- sampleSeed(1)

  # Create the command
  cmd = paste("java -jar", jar.path, as.character(msms.args), 
              "-ms", as.character(ms.args), "-seed", seed, ">", out.file)

  # Execute the command
  output <- system(cmd)

  if(!file.exists(out.file)) stop("msms simulation failed!")
  if(file.info(out.file)$size == 0) stop("msms output is empty!")
  
  return(out.file)
}

checkForMsms <- function(throw.error = TRUE, silent = FALSE) {
  if (isJaathaVariable('msms.jar')) {
    if (file.exists(getJaathaVariable('msms.jar'))) {  
      return(TRUE)
    }
  }

  # Works on Linux only maybe
  run.path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
  executables <- paste(c(run.path, getwd()), "/msms.jar", sep="")
  for (exe in executables) {
    if (file.exists(exe)) {
      if (!silent) message(paste("Using", exe, "as msms implementation\n"))
      setJaathaVariable('msms.jar', exe)     
      return(TRUE)
    }
  }

  if (throw.error) stop("No msms executable found.")
  return(FALSE)
}



# This function generates an string that contains an R command for generating
# an msms call to the current model.
generateMsmsOptionsCommand <- function(dm) {
  cmd <- c('c(')

  for (i in 1:dim(dm@features)[1] ) {
    type <- as.character(dm@features[i,"type"])
    feat <- unlist(dm@features[i, ])

    if (type == "pos.selection") {
      cmd <- c(cmd, '"-SI"', ',', feat['time.point'], ',', length(dm.getSampleSize(dm)), ',')
      start.freq <- rep(0, length(dm.getSampleSize(dm)))
      start.freq[ as.integer(feat['pop.source']) ] <- 0.0005
      cmd <- c(cmd, paste0('"', paste(start.freq, collapse=' '), '"'), ',')
      
      cmd <- c(cmd, '"-N 10000"', ',') 
      cmd <- c(cmd, '"-SA"', ',', feat['parameter'], ',') 
      cmd <- c(cmd, '"-Sp 0.5"', ',', '"-SForceKeep"', ',')
      cmd <- c(cmd, '"-threads 1"', ',')
    }
    
    else if (type == "bal.selection") {
      cmd <- c(cmd, '"-SI"', ',', feat['time.point'], ',', length(dm.getSampleSize(dm)), ',')
      start.freq <- rep(0, length(dm.getSampleSize(dm)))
      start.freq[ as.integer(feat['pop.source']) ] <- 0.0005
      cmd <- c(cmd, paste0('"', paste(start.freq, collapse=' '), '"'), ',')
      
      cmd <- c(cmd, '"-N 10000"', ',') 
      cmd <- c(cmd, '"-SAA 0"', ',',  '"-SAa"', ',', feat['parameter'], ',') 
      cmd <- c(cmd, '"-Sp 0.5"', ',', '"-SForceKeep"', ',')
      cmd <- c(cmd, '"-threads 1"', ',')
    }
  }

  cmd <- c(cmd, '" ")')
  cmd
}

createParameterEnv <- function(dm, parameters, ...) {
  par_env <- new.env()
  
  par.names <- dm.getParameters(dm)
  for (i in seq(along = par.names)){
    par_env[[ par.names[i] ]] <- parameters[i]
  }
  
  fixed.pars <- dm@parameters[dm@parameters$fixed, ]
  if (nrow(fixed.pars) > 0) {
    for (i in 1:nrow(fixed.pars)){
      par_env[[ fixed.pars$name[i] ]] <- fixed.pars$lower.range[i]
    }
  }
  
  additional_pars = list(...)
  for (i in seq(along = additional_pars)) {
    par_env[[names(additional_pars)[i]]] <- additional_pars[[i]]
  }
  
  par_env
}

generateMsmsOptions <- function(dm, parameters, locus) {
  msms.tmp <- createParameterEnv(dm, parameters, locus = locus)

  if ( !is.null( dm@options[['msms.cmd']] ) )
    cmd <- dm@options[['msms.cmd']]
  else
    cmd <- generateMsmsOptionsCommand(dm)
  cmd <- eval(parse(text=cmd), envir=msms.tmp)

  cmd
}

msmsSimFunc <- function(dm, parameters) {
  checkType(dm, "dm")
  checkType(parameters, "num")
  if (length(parameters) != dm.getNPar(dm)) 
    stop("Wrong number of parameters!")

  # Run all simulation in with one ms call if they loci are identical,
  # or call ms for each locus if there is variation between the loci.
  if (hasInterLocusVariation(dm)) {
    sim_reps <- 1:dm.getLociNumber(dm)
    sim_loci <- 1
  } else {
    sim_reps <- 1
    sim_loci <- dm.getLociNumber(dm)
  }
  
  # Run the simulation(s)
  msms.files <- lapply(sim_reps, function(locus) {
    ms.options <- paste(sum(dm.getSampleSize(dm)), sim_loci,
                        paste(generateMsOptions(dm, parameters, locus), 
                              collapse=" "))
    msms.options <- paste(generateMsmsOptions(dm, parameters, locus), 
                          collapse= " ")
    #print(c(ms.options, msms.options))
    callMsms(getJaathaVariable('msms.jar'), ms.options, msms.options)
  })
  
  # Parse the simulation output
  generateSumStats(msms.files, 'ms', parameters, dm)
}

finalizeMsms <- function(dm) {
  checkForMsms()
  dm@options[['ms.cmd']] <- generateMsOptionsCommand(dm)
  dm@options[['msms.cmd']] <- generateMsmsOptionsCommand(dm)
  return(dm)
}

printMsmsCommand <- function(dm) {
  msms.cmd <- printOptionsCmd(generateMsmsOptionsCommand(dm))
  ms.cmd <- printOptionsCmd(generateMsOptionsCommand(dm)) 
   
  cmd <- paste("msms", msms.cmd, 
               "-ms", sum(dm.getSampleSize(dm)), dm.getLociNumber(dm), ms.cmd)
  .print(cmd)
}

printOptionsCmd <- function(cmd) {
  cmd <- cmd[cmd != ","]
  cmd <- cmd[-c(1, length(cmd))]
  
  cmd <- paste(cmd, collapse=" ")
  
  cmd <- gsub(",", " ", cmd)
  cmd <- gsub('\"', "", cmd)
  cmd <- gsub('"', " ", cmd)  
}

#' @include dm_sim_program.R
createSimProgram("msms", possible.features, possible.sum.stats,
                 msmsSimFunc, finalizeMsms, printMsmsCommand, priority=40)
