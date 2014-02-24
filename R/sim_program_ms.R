if (!exists('sim.progs')) sim.progs <- list()
.ms <- new.env(parent=emptyenv())
sim.progs$ms <- .ms

# Supported features
.ms$features <- c("mutation", "migration", "split",
                  "recombination", "size.change", "growth")

# Supported summary statistics
.ms$sum.stats <- c("jsfs", "4pc", "tree", "seg.sites")

# Function to perform the actual simulations with ms
.ms$simulate <- function(dm, parameters) {
  checkType(dm, "dm")
  checkType(parameters, "num")

  if (length(parameters) != dm.getNPar(dm)) stop("Wrong number of parameters!")

  ms.options <- generateMsOptions(dm, parameters)
  ms.out <- callMs(ms.options, dm)

  sum.stats <- list(pars=parameters)

  if ("jsfs" %in% dm@sum.stats) {
    sum.stats[['jsfs']] <- msOut2Jsfs(dm, ms.out)
  }

  unlink(ms.out)
  return(sum.stats)
}

.ms$finalize <- function(dm) {
  dm@options[['ms.cmd']] <- generateMsOptionsCommand(dm)
  return(dm)
}

## Function to perform simulation using ms 
## 
## @param opts The options to pass to ms. Must either be a character or character
## vector.
## @param dm The demographic model we are using
## @return The file containing the output of ms
.ms$callMs <- function(opts, dm) {
  if (missing(opts)) stop("No options given!")
  opts <- unlist(strsplit(opts, " "))

  ms.file <- getTempFile("ms")

  ms(sum(dm@sampleSizes), dm@nLoci, opts, ms.file)
  return(ms.file)
}


# This function generates an string that contains an R command for generating
# an ms call to the current model.
.ms$generateMsOptionsCommand <- function(dm) {
  nSample <- dm@sampleSizes
  cmd <- c('c(')
  cmd <- c(cmd,'"-I"', ",", length(nSample), ',', 
           paste(nSample, collapse=","), ',')

  for (i in 1:dim(dm@features)[1] ) {
    type <- as.character(dm@features[i,"type"])
    feat <- unlist(dm@features[i, ])

    if (type == "mutation") {
      cmd <- c(cmd,'"-t"', ',', feat["parameter"], ',')
    }

    if (type == "split") {
      cmd <- c(cmd, '"-ej"', ',', feat["time.point"], ',',
               feat["pop.sink"], ',', feat["pop.source"], ',')
    }

    if (type == "migration")
      cmd <- c(cmd, '"-em"', ',', feat['time.point'], ',',
               feat['pop.sink'], ',', feat['pop.source']  , ',',
               feat['parameter'], ',')

    if (type == "recombination") 
      cmd <- c(cmd, '"-r"', ',', feat['parameter'], ',', dm@seqLength, ',')

    if (type == "size.change"){
      cmd <- c(cmd, '"-en"', ',', feat['time.point'], ',',
               feat["pop.source"], ',', feat['parameter'], ',')
    }

    if (type == "growth"){
      cmd <- c(cmd, '"-eg"', ',' , feat["time.point"], ',',
               feat["pop.source"], ',', feat["parameter"], ',')
    }
  }

  cmd <- c(cmd, '"-T")')
}

.ms$generateMsOptions <- function(dm, parameters) {
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

.ms$printMsCommand <- function(dm) {
  cmd <- generateMsOptionsCommand(dm)

  cmd <- cmd[cmd != ","]
  cmd <- cmd[-c(1, length(cmd))]

  cmd <- paste(cmd, collapse=" ")

  cmd <- gsub(",", " ", cmd)
  cmd <- gsub('\"', "", cmd)
  cmd <- gsub('"', " ", cmd)

  return(cmd)
}

.ms$msOut2Jsfs <- function(dm, ms.out) {
  jsfs <- matrix(.Call("msFile2jsfs", ms.out, dm@sampleSizes[1], 
                       dm@sampleSizes[2]),
                 dm@sampleSizes[1] + 1 ,
                 dm@sampleSizes[2] + 1,
                 byrow=T)
  return(jsfs)
}
