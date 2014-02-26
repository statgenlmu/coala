#---------------------------------------------------------------
# DemographicModel.R
# Class for representing a model of the evolutionary development
# of two different species.
# 
# Authors:  Paul R. Staab & Lisha Mathew 
# Email:    staab ( at ) bio.lmu.de
# Date:     2013-09-04
# Licence:  GPLv3 or later
#--------------------------------------------------------------

#' @docType package
#' @name EvoModel-package
#' @title Fast estimation of demographic parameters
#' @keywords package
#' @importFrom phyclust ms
#' @import methods
#' @import Rcpp
#' @useDynLib EvoModel
NULL

init <- function() {
  features <<- data.frame(type=character(),
                                 parameter=character(),
                                 pop.source=numeric(),
                                 pop.sink=numeric(),
                                 time.point=character(),
                                 group=numeric(),
                                 stringsAsFactors=F )

  parameters <<- data.frame(parameter=character(),
                                   fixed=logical(),
                                   lower.range=numeric(),
                                   upper.range=numeric(),
                                   stringsAsFactors=F )

  sum.stats <<- c("jsfs")
  options <<- list()
}


#' The EvoModel Class
#' 
#' This class represents a model of the evolution of biological sequences. 
#'
#' @section Slots:
#' \describe{
#'   \item{\code{features}:}{A data.frame that contains the features the model has.}
#'   \item{\code{parameters}:}{A data.frame that contains the model's parameters.}
#'   \item{\code{sum.stats}:}{A vector with the names of Summary Statistics that
#'                            will be created if the model is simulated.}
#'   \item{\code{sim.prog}:}{The program that is used for simulation.}
#'   \item{\code{options}:}{A list of additional options that the program may
#'                          use.}
#' }
#' @name EvoModel
#' @rdname EvoModel
#' @aliases EvoModel-class
#' @author Paul R. Staab
#' @export
EvoModel <- setRefClass("EvoModel",
         fields = list(features="data.frame",
                       parameters="data.frame",
                       sum.stats="character",
                       sim.prog="environment",
                       options="list"))

rm(init)

#-----------------------------------------------------------------------
# Initialization
#-----------------------------------------------------------------------





#-----------------------------------------------------------------------
# Print
#-----------------------------------------------------------------------

.show <- function(object){
  dm <- object

  features <- dm@features[!is.na(dm@features$parameter),]

  if (nrow(features) == 0) {
    cat("Your demographic model has no features so far.\n")
    return()
  }

  cat("Your demographic model has the following parameters:\n")
  for (rw in 1:nrow(features)) {
    type <- features[rw,'type'] 
    lR <- as.character(features[rw,'lower.range']) 
    uR <- as.character(features[rw,'upper.range'])
    parname <- features[rw,'parameter']

    if (type == "split")
      cat("-",parname,": A split into two pop.sources between",lR,"and",uR,
          "Ne generations ago.\n")
    else if (type == "mutation")
      cat("-",parname,": A scaled mutation rate between",lR,"and",uR,"\n")
    else if (type == "recombination")
      cat("-",parname,": A scaled recombination rate between",lR,"and",uR,"\n")
    else if (type == "migration")
      cat("-",parname,":  A scaled migration rate between",lR,"and",uR,"\n")
    else if (type == "presentSize")
      cat("-",parname,":  At sample time, the second pop.source two was between",
          lR,"and",uR,"times as the first one\n")
    else if (type == "splitSize") 
      cat("-",parname,":  At time of the pop.source split, the second",
          "pop.source two was between",lR,"and",uR,"times as the first one\n")
    else if (type == "migration")
      cat("-",parname,":  A scaled migration rate between",lR,"and",uR,"\n")
    else {
      cat("- Other:\n")
      print(dm@features[rw,])
    }
  }


  features <- dm@features[is.na(dm@features$parameter),]
  if (dim(features)[1] > 0){
    cat("\n")
    cat("Additionally, it has the following features:\n")
    for (rw in 1:(dim(features)[1])) {
      type <- features[rw,'type']
      lR <- as.character(features[rw,'lower.range']) 

      if (type == "split")
        cat("- A split into two pop.sources",lR,"Ne generations ago.\n")
      else if (type == "mutation")
        cat("- A fixed scaled mutation rate of",lR,"\n")
      else if (type == "recombination")
        cat("- A fixed scaled recombination rate of",lR,"\n")
      else if (type == "migration")
        cat("- A fixed scaled migration rate of",lR,"\n")
      else if (type == "presentSize")
        cat("- At sample time, the second pop.source two was",
            lR,"times as the first one\n")
      else if (type == "splitSize") 
        cat("- At time of the pop.source split, the second pop.source
            two was",lR,"times as the first one\n")
          else {
            cat("- Other:\n")
            print(dm@features[rw,])
          }
    }
  }
} 
#setMethod("show","DemographicModel",.show)
rm(.show)



#------------------------------------------------------------------------------
# Private functions
#------------------------------------------------------------------------------

#---------------------------------------------------------------------
# addParameter()
#---------------------------------------------------------------------
#' Create a parameter that can be used for one or more features
#'
#' The function creates a new model parameter. It can either have a fixed value
#' or you can enter a range of possible values if you want to estimate it.
#'
#' @param dm    The demographic model to which the parameter will be added
#' @param par.name The name of the parameter. You can use this name later to
#'                 access the parameter from a model feature
#' @param lower.boundary The lower boundary of the range within which the
#'                       parameter value will be estimated. Don't specify
#'                       'fixed.value' if you want to do so.
#' @param upper.boundary Like 'lower.boundary', but the upper end of the
#'                       parameter range.
#' @param fixed.value    If this argument is given, than rather than being estimated 
#'                       a fixed value will be used.
#' @return The original model extended with the new parameter.
#' @export
#' @examples
#' dm <- dm.createThetaTauModel(c(15,23), 100)
#' dm <- dm.addParameter(dm, "mig", 0.1, 5)
#' dm <- dm.addMigration(dm, par.new=FALSE, parameter="mig", pop.from=1, pop.to=2)
#' dm <- dm.addMigration(dm, par.new=FALSE, parameter="2*mig", pop.from=2, pop.to=1)
dm.addParameter <- function(dm, par.name, lower.boundary, upper.boundary, fixed.value) {

  if (missing(lower.boundary)) lower.boundary <- NA
  if (missing(upper.boundary)) upper.boundary <- NA
  if (missing(fixed.value)) fixed.value <- NA

  checkType(dm,          c("dm",   "s"), T, F)
  checkType(par.name,    c("char", "s"), T, F)
  checkType(lower.boundary, c("num",  "s"), T, T)
  checkType(upper.boundary, c("num",  "s"), T, T)
  checkType(fixed.value, c("num",  "s"), T, T)

  if (par.name %in% dm.getParameters(dm, T))
    stop("There is already a parameter with name ",par.name)

  if ( !is.na(fixed.value) ) {
    lower.boundary <- fixed.value
    upper.boundary <- fixed.value
  }

  fixed <- F
  if (is.na(lower.boundary) | is.na(upper.boundary)) {
    fixed <- T
  } else {
    if (lower.boundary == upper.boundary) {
      fixed <- T
    }
  }
 
  dm <- appendToParameters(dm, par.name, fixed, lower.boundary, upper.boundary)

  dm <- makeThetaLast(dm)
  return(dm)
}

#' Low level function for adding a new feature to a demographic Model
#'
#' Use this function to add a feature to a dm if there is no higher level
#' "dm.add*"-function availible.
#'
#' @param dm       The demographic model to which the feature will be added 
#' @param type     The type of the feature coded as a character
#' @param lower.range The lower boundary for the value of the parameter
#' @param upper.range The upper boundary for the value of the parameter
#' @param fixed.value If given, the parameter will be set to a fixed value. This
#'                 is equivalent to seting lower.range equal to upper.range.
#' @param pop.source The source population if availible (think e.g. of migration)
#' @param pop.sink   The target or "sink" population if availible (think e.g. of migration)
#' @param time.point Normally the point in backwards time where the feature
#'                   starts.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'lower.range' and 'upper.range' or 
#'            'fixed.value'. It will be named 'parameter'.
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param parameter Either the name of the parameter (par.new=TRUE), or an R expression
#'            possibly containing one or more previously created parameter
#'            names.
#' @param group     For genomic features, different groups can be created.
#' @return          The extended demographic model.
addFeature <- function(dm, type, parameter=NA,
                       lower.range=NA, upper.range=NA, fixed.value=NA, par.new=T,
                       pop.source=NA, pop.sink=NA,
                       time.point=NA, group=NA) {
  
  if (missing(par.new))     par.new <- T
  if (missing(parameter))   parameter <- NA
  if (missing(pop.source))  pop.source <- NA
  if (missing(pop.sink))    pop.sink <- NA
  if (missing(time.point))  time.point <- NA
  if (missing(group))       group <- NA

  # Check inputs
  checkType(dm,          c("dm",   "s"), T, F)
  checkType(type,        c("char", "s"), T, F)
  checkType(par.new,     c("bool", "s"), T, F)
  checkType(parameter,   c("char", "s"), T, T)
  checkType(pop.source,  c("num",  "s"), T, T)
  checkType(pop.sink,    c("num",  "s"), T, T)
  checkType(time.point,  c("char", "s"), T, T)
  checkType(group,       c("num",  "s"), T, T)

  if (par.new) dm <- dm.addParameter(dm, parameter, lower.range, upper.range, fixed.value)

  # Append the feature
  dm <- appendToFeatures(dm = dm,
                         type = type,
                         parameter = parameter,
                         pop.source = pop.source,
                         pop.sink = pop.sink,
                         time.point = time.point,
                         group = group)

  # Update some technical properties of the model
  dm <- .dm.selectSimProg(dm)
}

#' @title Add Summary Statistic
#' 
#' @param sum.stat The name of a summary statistic.
#' 
#' @name EvoModel_addSummaryStatistic
NULL
EvoModel$methods(addSummaryStatistic = function(sum.stat) {
  'see \\code{\\link{EvoModel_addSummaryStatistic}}.'
  checkType(dm, "dm")
  checkType(sum.stat, "char")
  dm@sum.stats <- c(dm@sum.stats, sum.stat)
  dm <- .dm.selectSimProg(dm)
  return(dm)
})


# Helper function that appends a feature to the "feature" dataframe
# Does not check the feature for consistency
# This should only be used by addFeature().
appendToFeatures <- function(dm, type, parameter, pop.source, 
                             pop.sink, time.point, group    ) {
  
  new.feature <- data.frame(type=type,
                            parameter=parameter,
                            pop.source=pop.source,
                            pop.sink=pop.sink,
                            time.point=time.point,
                            group=group,
                            stringsAsFactors=F)

  dm@features <- rbind(dm@features, new.feature)
  return(dm)
}

# Helper function that appends a parameter to the "parameters" dataframe
# Does not check the feature for consistency
# This should only be used by addFeature().
appendToParameters <- function(dm, name, fixed, lower.range, upper.range) {

  new.parameter <- data.frame(name=name,
                              fixed=fixed,
                              lower.range=lower.range,
                              upper.range=upper.range,
                              stringsAsFactors=F)
  
  dm@parameters <- rbind(dm@parameters, new.parameter)
  return(dm)
}


# Gets the availible populations
getPopulations <- function(dm){
  # Every population other than the ancestral (=0) should be listed
  # the pop.sink field of its speciation event.
  populations <- c(1, dm@features$pop.sink)
  populations <- unique(populations[!is.na(populations)])
}


# To use Watersons estimator, Jaatha requires that the mutation parameter
# is the last parameter. This swishes it to the end.
makeThetaLast <- function(dm) {
  checkType(dm, c("dm", "s"), T, F)
  
  par.name <- dm@features$parameter[dm@features$type == "mutation"]
  if (length(par.name) == 0) return(dm)
  if (length(par.name) > 1) stop("Multiple Mutations features found")

  pars <- dm@parameters
  mut.line <- (1:nrow(pars))[pars$name == par.name]
  last.line <- nrow(pars)
  if (mut.line == last.line) return(dm)

  new.order <- 1:last.line
  new.order[c(mut.line, last.line)] <- c(last.line, mut.line)

  dm@parameters <- pars[new.order, ]
  return(dm)
}


# Checks if a vector of parameters is within the ranges of the model
checkParInRange <- function(dm, param) {
  if (length(param) != dm.getNPar(dm)) stop("Wrong number of parameters")

  ranges <- dm.getParRanges(dm)
  in.range <- all(ranges[, 1]-1e-11 <= param & param <= ranges[, 2]+1e-11)
  if (!in.range) stop("Parameter combination out of range")
}

# Selects a program for simulation that is capable of all current features
.dm.selectSimProg <- function(dm) {
  for (i in seq(along = .jaatha$simProgs)){
    if (all(dm@features$type %in% .jaatha$simProgs[[i]]@possible.features) & 
        all(dm@sum.stats %in% .jaatha$simProgs[[i]]@possible.sum.stats)) {
      dm@currentSimProg <- .jaatha$simProgs[[i]]
      .log2("Using", dm@currentSimProg@name, "for simulations")
      return(dm)
    }
  }
  stop("No suitable simulation software found!")
  return(dm)
}

finalizeDM <- function(dm) {
  simProg <- dm@currentSimProg
  .log2("Finalizing demographic model. Simprog:", simProg@name)
  dm <- simProg@finalizationFunc(dm)
  return(dm)
}


#------------------------------------------------------------------------------
# Getters & Setters for Jaatha
#------------------------------------------------------------------------------
dm.getParameters <- function(dm, fixed=F) {
  checkType(dm,"dm")
  param <- dm@parameters$name
  param <- param[dm@parameters$fixed == fixed]
  return(param)
}

dm.getNPar <- function(dm){
  checkType(dm, "dm")
  return(length(dm.getParameters(dm)))
}

dm.getParRanges <- function(dm, inklExtTheta=T){
  #checkType(dm,"dm")
  par.ranges <- dm@parameters[!dm@parameters$fixed, c("lower.range","upper.range")]
  rownames(par.ranges) <- dm@parameters[!dm@parameters$fixed,"name"]
  if (!inklExtTheta) par.ranges <- par.ranges[1:dm.getNPar(dm),]
  return(par.ranges)
}

fixTheta <- function(dm) {
  param <- dm@parameters
  param[nrow(param), "fixed"] <- T
  dm@parameters <- param
  return(dm)
}

getThetaName <- function(dm){
  return(dm@parameters[nrow(dm@parameters), 'name'])
}

getThetaRange <- function(dm){
  return(dm@parameters[nrow(dm@parameters),c('lower.range', 'upper.range')])
}





#------------------------------------------------------------------------------
# Creation new models
#------------------------------------------------------------------------------

#---------------------------------------------------------------------
# dm.createDemographicModel()
#---------------------------------------------------------------------
#' Create a basic demographic model
#' 
#' This function creates a basic empty demographic model, which
#' is returned. Features like mutation, pop.source splits and 
#' migration can be added afterwards.
#'
#' @param sample.sizes Number of individuals that are sampled. If your model 
#'            consists of multiple populations, this needs to be a vector
#'            containing the sample sizes from each population.
#' @param loci.num     Number of loci that will be simulated
#' @param seq.length   (Average) number of bases for each locus
# @param finiteSites If 'TRUE', a finite sites mutation model is assumed
#            instead of an infinite sites one.
# @param tsTvRatio   Transition transversion ratio
#' @param log.level   An integer specifing the amount of user readable output 
#'                    that will be produced.
#'                    0 = no output, 1 = normal verbosity, ..., 
#'                    3 = maximal verbosity.
#' @param log.file    If set, the debug output will be written into the given file
#' @return            The demographic model
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sample.sizes=c(25,25), loci.num=100)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.createDemographicModel <- function(sample.sizes, loci.num, seq.length=1000, 
                                      #finiteSites=F, tsTvRatio=.33, 
                                      log.level, log.file) {
  setLogging(log.level, log.file)
  dm <- new("DemographicModel", sample.sizes, loci.num, seq.length, F, .33)
  return(dm)
}



#---------------------------------------------------------------------------------
# Front end functions for adding features
#---------------------------------------------------------------------------------

#--------------------------------------------------------------------
# dm.addMutation()
#--------------------------------------------------------------------
#' Adds mutations to a demographic model
#' 
#' This functions adds the assumption to the model that neutral mutations
#' occur in the genomes at a constant rate. The rate is quantified through
#' a parameter usually named theta in population genetics. It equals 4*Ne*mu,
#' where Ne is the (effective) number of diploid individuals in the ancestral
#' population and mu is the neutral mutation rate for an entire locus.
#'
#' @param dm  The demographic model to which mutations should be added
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'lower.range' and 'upper.range' or 
#'            'fixed.value'. It will be named 'new.par.name'.
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param lower.range  If you want to estimate the mutation rate, this 
#'            will be used as the smallest possible value.
#' @param upper.range  If you want to estimate the mutation rate, this 
#'            will be used as the largest possible value.
#' @param fixed.value  If specified, the mutation rate will not be 
#'            estimated, but assumed to be fixed at the given value.
#' @param new.par.name The name for the new parameter.
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "theta" will use
#'            an parameter with name theta that you have previously 
#'            created. You can also use R expression here, so "2*theta"
#'            or "5*theta+2*tau" (if tau is another parameter) will also
#'            work (also it does not make much sense).
#' @return    The demographic model with mutation.
#' @export
#'
#' @examples
#' # Create a new parameter
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addSpeciationEvent(dm, 0.01, 5)
#' dm <- dm.addMutation(dm, 1, 20)
#'
#' # Create a new fixed parameter
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addSpeciationEvent(dm, 0.01,5)
#' dm <- dm.addMutation(dm, fixed.value=7)
#'
#' # Use an existing parameter
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addParameter(dm, "theta", 0.01, 5)
#' dm <- dm.addMutation(dm, par.new=FALSE, parameter="2*log(theta)+1")
dm.addMutation <- function(dm, lower.range, upper.range, fixed.value,
                           par.new=T, new.par.name="theta", parameter) {

  if ( missing(lower.range) & missing(upper.range) & 
       missing(fixed.value) & missing(parameter) ) {
    stop("Either a parameter range or a fixed value is required. If you 
          want to use 'externalTheta', please use Version 2.0.2")
  }

  if (par.new) parameter <- new.par.name
  dm <- addFeature(dm, "mutation", parameter, lower.range, upper.range,
                   fixed.value, par.new=par.new, time.point=NA)
  return(dm)
}


#-------------------------------------------------------------------
#  dm.addRecombination()
#-------------------------------------------------------------------
#' Adds recombination events to a demographic model
#'
#' This function add the assumption to the model that recombination
#' events may occur within each locus. The corresponding parameter
#' - usually name rho - equals 4*Ne*r, where Ne is the number of
#' diploid individuals in the ancestral population and r is the 
#' probability that a recombination event within the locus will
#' occur in one generation. Even when using an infinite sites
#' mutation model, this assumes an finite locus length which is given
#' by the 'seq.length' parameter of the demographic model.
#'
#' Please note that it does not make sense to estimate recombination
#' rates with Jaatha because it assumes unlinked loci.
#' 
#' @param dm  The demographic model to which recombination events should be added.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'lower.range' and 'upper.range' or 
#'            'fixed.value'. It will be named 'new.par.name'.
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param lower.range  If you want to estimate the recombination rate (see note
#             above), this will be used as the smallest possible value.
#' @param upper.range  Same as lower.range, but the largest possible value.
#' @param fixed.value  If specified, the mutation rate will not be estimated,
#'                     but assumed to have the given value.
#' @param new.par.name The name for the new parameter.
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "rho" will use
#'            an parameter with name theta that you have previously 
#'            created. You can also use R expression here, so "2*rho"
#'            or "5*rho+2*tau" (if tau is another parameter) will also
#'            work (also it does not make much sense).
#' @return    The demographic model with recombination
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addSpeciationEvent(dm, 0.01, 5)
#' dm <- dm.addRecombination(dm, fixed=20)
#' dm <- dm.addMutation(dm, 1, 20)
dm.addRecombination <- function(dm, lower.range, upper.range, fixed.value,
                                par.new=T, new.par.name="rho", parameter) {

  if (par.new) parameter <- new.par.name
  dm <- addFeature(dm, "recombination", parameter, lower.range, upper.range,
                   fixed.value, par.new=par.new, time.point="0")
  return(dm)
}


#-------------------------------------------------------------------
#  dm.addMigration
#-------------------------------------------------------------------
#' Add migration/gene flow between two populations to a demographic model
#'
#' This function adds the assumption to the model that some individuals
#' 'migrate' from one sub-population to another, i.e. they leave the one
#' and become a member of the other. This is usually used to model ongoing
#' gene flow through hybridisation after the populations separated. 
#'
#' You can enter a time ('time.start') at which the migration is 
#' assumed to start (looking backwards in time). From that time on, a 
#' fixed number of migrants move from population 'pop.from' to
#' population 'pop.to' each generation. This number is given via this 
#' feature's parameter, which equals 4*Ne*m,  where Ne is the number of
#' diploid individuals in the ancestral population and m is the 
#' fraction of 'pop.to' that is replaced with migrants each generation. 
#' If 'pop.to' has also size Ne, than this is just the
#' expected number of individuals that migrate each generation.
#' 
#' You can add different mutation rates at different times to your model.
#' Then each rate will be used for the period from its time point to
#' the next. Migration from and to an population always ends with the 
#' speciation event in which the population is created.
#'
#' @param dm  The demographic model to which recombination events should be added.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'lower.range' and 'upper.range' or 
#'            'fixed.value'. It will be named 'new.par.name'.
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param lower.range  If you want to estimate the migration parameter (see note
#             above), this will be used as the smallest possible value.
#' @param upper.range  Same as lower.range, but the largest possible value.
#' @param fixed.value  If specified, the migration rate will not be estimated,
#'                     but assumed to have the given value.
#' @param new.par.name The name for the new parameter.
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "M" will use
#'            an parameter with name M that you have previously 
#'            created. You can also use R expression here, so "2*M"
#'            or "5*M+2*tau" (if tau is another parameter) will also
#'            work (also this does not make much sense).
#' @param pop.from The population from which the individuals leave.
#' @param pop.to The population to which the individuals move.
#' @param time.start The time point at which the migration with this rate
#'            starts.
#' @return    The demographic model with migration
#' @export
#'
#' @examples
#' # Constant asymmetric migration
#' dm <- dm.createThetaTauModel(c(25,25), 100)
#' dm <- dm.addMigration(dm, 0.01, 5, pop.from=1, pop.to=2, time.start="0")
#'
#' # Stepwise decreasing mutation
#' dm <- dm.createThetaTauModel(c(25,25), 100)
#' dm <- dm.addMigration(dm, 0.01, 5, pop.from=1, pop.to=2, new.par.name="M", 
#'                       time.start="0")
#' dm <- dm.addMigration(dm, pop.from=1, pop.to=2, par.new=FALSE, 
#'                       parameter="0.5*M", time.start="0.5*tau")
dm.addMigration <- function(dm, lower.range, upper.range, fixed.value,
                            pop.from, pop.to, time.start="0",
                            par.new=T, new.par.name="M",
                            parameter) {

  checkType(pop.from, "num", T, F)
  checkType(pop.to,   "num", T, F)

  if (par.new) parameter <- new.par.name
  dm <- addFeature(dm, "migration", parameter, lower.range, upper.range,
                   fixed.value, par.new, pop.from, pop.to, time.start)
  return(dm)
}


#-------------------------------------------------------------------
#  dm.addSymmetricMigration
#-------------------------------------------------------------------
#' Adds symmetric migration between all populations
#'
#' This adds migration between all subpopulation, all with the same
#' rate. Please look at the documentation for \link{dm.addMigration} for 
#' detailed information about migration.
#' 
#' @param dm  The demographic model to which migration events should be added.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'lower.range' and 'upper.range' or 
#'            'fixed.value'. It will be named 'new.par.name'.
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param lower.range  If you want to estimate the migration parameter (see note
#             above), this will be used as the smallest possible value.
#' @param upper.range  Same as lower.range, but the largest possible value.
#' @param fixed.value  If specified, the migration rate will not be estimated,
#'                     but assumed to have the given value.
#' @param new.par.name The name for the new parameter.
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "M" will use
#'            an parameter with name M that you have previously 
#'            created. You can also use R expression here, so "2*M"
#'            or "5*M+2*tau" (if tau is another parameter) will also
#'            work (also this does not make much sense).
#' @param time.start The time point at which the migration with this rate
#'            starts.
#' @return    The demographic model with migration
#' @export
#'
#' @examples
#' dm <- dm.createThetaTauModel(c(25,25), 100)
#' dm <- dm.addSymmetricMigration(dm, 0.01, 5, time.start="0.5*tau")
dm.addSymmetricMigration <- function(dm, lower.range, upper.range, fixed.value,
                                     par.new=T, new.par.name="M", parameter, 
                                     time.start="0") {

  if (par.new) dm <- dm.addParameter(dm, new.par.name, lower.range, upper.range, fixed.value)
  if (par.new) parameter <- new.par.name

  for (i in getPopulations(dm)) {
    for (j in getPopulations(dm)) {
      if (i==j) next
      dm <- dm.addMigration(dm, par.new=F, parameter=parameter, time.start=time.start,
                            pop.from=i, pop.to=j)
    }
  }

  return(dm)
}


#-------------------------------------------------------------------
# dm.addSpeciationEvent
#-------------------------------------------------------------------
#' Adds a speciation event to a demographic model
#'
#' You can use this function the create a new population that splits 
#' of from an existing population at a given time in the past. The time
#' can be given as parameter or as an expression based on previously
#' generated time points.
#'
#' As always, time in measured in Number of 4Ne generations in the past,
#' where Ne is the (effective) size of the ancestral population.
#'
#' The command will print the number of the new population, which will
#' be the number of previously existing populations plus one. The ancestral
#' population has number "1".
#'
#' @param dm  The demographic model to which the split should be added.
#' @param new.time.point  If 'TRUE' a new parameter will be created using the
#'            arguments 'min.time' and 'max.time' or 
#'            'fixed.time'. It will be named 'new.time.point.name'.
#'            If 'FALSE' the argument 'time.point' 
#'            will be evaluated instead.
#' @param min.time  If you want to estimate the time point, this will be 
#'            used as the smallest possible value.
#' @param max.time  Same as min.time, but the largest possible value.
#' @param fixed.time  If specified, the time.point will not be estimated,
#'            but assumed to have the given value.
#' @param new.time.point.name The name for the new time point.
#' @param in.population The number of the population in which the spilt
#'            occurs. See above for more information.
#' @param time.point  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "tau" will use
#'            an parameter with name tau that you have previously 
#'            created. You can also use R expression here, i.e. "2*tau"
#'            or "5*M+2*tau" (if M is another parameter) will also
#'            work (also this does not make much sense).
#' @return    The demographic model with a split.
#' @export
#' @examples
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
dm.addSpeciationEvent <- function(dm, min.time, max.time, fixed.time, 
                                  in.population=1, 
                                  new.time.point=T, new.time.point.name=NA,
                                  time.point=NA){

  checkType(dm,                  c("dm",  "s"),  T, F)
  checkType(in.population,       c("num", "s"),  T, F)
  checkType(min.time,            c("num", "s"),  F, F)
  checkType(max.time,            c("num", "s"),  F, F)
  checkType(fixed.time,          c("num", "s"),  F, F)
  checkType(new.time.point.name, c("char","s"),  F, T)
  checkType(time.point,          c("char", "s"),  F, T)
  
  # in.population valid?
  populations <- getPopulations(dm)
  if (!in.population %in% populations) 
    stop("Population", in.population, "not known")

  # time.points
  if (new.time.point) {
    if(is.na(new.time.point.name)) 
      new.time.point.name <- paste("t_split_", in.population, sep="")
    dm <- dm.addParameter(dm, new.time.point.name, min.time, 
                          max.time, fixed.time)
    time.point <- new.time.point.name
  } else {
    # Check if valid
  }

  new.pop <- max(populations) + 1

  return(addFeature(dm, "split", par.new = F, 
                    pop.source=in.population,
                    pop.sink=new.pop,
                    time.point=time.point))
}


#-------------------------------------------------------------------
# dm.addSizeChange
#-------------------------------------------------------------------
#' Adds an instantaneous change of the population size of one 
#' population to a model.
#'
#' This function changes the effective population size of one 
#' population. The change is performed at a given time point
#' ('at.time') and applies to the time interval farther into 
#' the past from this point. The population size is set to a
#' factor of the size of the ancestral population Ne.
#'
#' If you want to add a slow, continuous change over some time,
#' then use the \link{dm.addGrowth} function.
#' 
#' @param dm  The demographic model to which the size change should be added.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'min.size.factor' and 'max.size.factor' or 
#'            'fixed.size.factor'. It will be named 'new.time.point.name'
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param min.size.factor  If you want to estimate the size factor, this will be 
#'            used as the smallest possible value.
#' @param max.size.factor  Same as min.size.factor, but the largest possible value.
#' @param fixed.size.factor  If specified, the size factor not be estimated,
#'            but assumed to have the given value.
#' @param new.par.name Name for the new parameter.
#' @param population The number of the population in which the spilt
#'            occurs. See \link{dm.addSpeciationEvent} for more information.
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "tau" will use
#'            an parameter with name tau that you have previously 
#'            created. You can also use R expression here, i.e. "2*tau"
#'            or "5*M+2*tau" (if M is another parameter) will also
#'            work (also this does not make much sense).
#' @param at.time The time point at which the size changes.
#' @return    The demographic model with a size change.
#' @export
#' @examples
#' # A model with one smaller population
#' dm <- dm.createThetaTauModel(c(20,37), 88)
#' dm <- dm.addSizeChange(dm, 0.1, 1, population=2, at.time="0")
dm.addSizeChange <- function(dm, min.size.factor, max.size.factor,
                             fixed.size.factor, par.new=T, new.par.name="q",
                             parameter, 
                             population, at.time="0") {

  checkType(population, c("num",  "s"), T, F)
  checkType(at.time,    c("char", "s"), T, F)

  if (par.new) parameter <- new.par.name

  dm <- addFeature(dm, "size.change", parameter, min.size.factor, max.size.factor,
                   fixed.size.factor, par.new, population, NA, at.time)

  return(dm)
}


#-------------------------------------------------------------------
# dm.addGrowth
#-------------------------------------------------------------------
#' Adds an growth or decline of the population size of one 
#' population to a model.
#'
#' This function changes the growth factor of a population at given 
#' point in time ('at.time'). This factor than applies to the time 
#' interval farther into the past from this point. 
#'
#' The population size changes by factor exp(-alpha*t), where alpha
#' is the growth parameter and t is the time since the growth has 
#' started. Hence, for positive alpha, the population will 'decline 
#' backwards in time' or grow forwards in time. Similar, will decline
#' in forwards time for a negative value of alpha.
#'
#' If you want to add an instantaneous change of the population size,
#' then use the \link{dm.addSizeChange} function.
#' 
#' @param dm  The demographic model to which the size change should be added.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'min.growth.rate' and 'max.growth.rate' or 
#'            'fixed.growth.rate'. It will be named 'new.par.name'
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param min.growth.rate  If you want to estimate the growth rate, this will be 
#'            used as the smallest possible value.
#' @param max.growth.rate  Same as min.growth.rate, but the largest possible value.
#' @param fixed.growth.rate  If specified, the growth rate will not be estimated,
#'            but assumed to have the given value.
#' @param new.par.name Name for the new parameter.
#' @param population The number of the population in which the spilt
#'            occurs. See \link{dm.addSpeciationEvent} for more information.
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "alpha" will use
#'            an parameter with name tau that you have previously 
#'            created. You can also use R expression here, i.e. "2*alpha"
#'            or "5*M+2*alpha" (if M is another parameter) will also
#'            work (also the latter does not make much sense).
#' @param at.time The time point at which the size changes.
#' @return    The demographic model with a size change.
#' @export
#' @examples
#' # A model with one smaller population
#' dm <- dm.createThetaTauModel(c(20,37), 88)
#' dm <- dm.addGrowth(dm, 0.1, 2, population=2, at.time="0")
dm.addGrowth <- function(dm, min.growth.rate, max.growth.rate, fixed.growth.rate, 
                         par.new=T, new.par.name="alpha", parameter, 
                         population, at.time="0") {

  checkType(population, c("num",  "s"), T, F)
  checkType(at.time,    c("char", "s"), T, F)

  if (par.new) parameter <- new.par.name

  dm <- addFeature(dm, "growth", parameter, min.growth.rate, max.growth.rate,
                   fixed.growth.rate, par.new, population, NA, at.time)

  return(dm)
}


#-------------------------------------------------------------------
# dm.setMutationModel
#-------------------------------------------------------------------
#' Defines what mutation model is used for simulations
#'
#' As default, we simulate mutation using the Infinite Sites Model. 
#' Using the function, you can change it either to the Hasegawa, Kishino and
#' Yano (HKY), to the Felsenstein and Churchill 96 (F84) or to the Generalised 
#' time reversible (GTR) model. This requires that seq-gen is installed on our system.
#'
#' The HKY and F84 models use the the arguments 'base.frequencies' and
#' 'tstv.ratio'. The GTR model uses 'gtr.rates'.
#'
#' @param dm  The demographic model for which the mutation model will be set.
#' @param mutation.model  The mutation model you want to use. Can be HKY, F84 or GTR.
#' @param tstv.ratio The ratio of transitions to transversions. The default is
#'                   0.5, which means that all amino acid substitutions are
#'                   equally likely. In this case, the HKY model is identical to
#'                   the Felsenstein 81 model.
#' @param base.frequencies The equilibrium frequencies of the four bases. 
#'                   Must be a numeric vector of length four. 
#'                   Order is A, C, G, T.
#' @param gtr.rates  The rates for the amino acid substitutions. Must be a
#'                   numeric vector of length six. Order: A->C, A->G, A->T, C->G, C->T, G->T.
#' @return    The demographic model with the new mutation model.
#' @export
#' @examples
#' dm <- dm.createThetaTauModel(10:11, 10, 100)
#' dm <- dm.addOutgroup(dm, "2*tau")
#' dm.hky <- dm.setMutationModel(dm, "HKY", c(0.2, 0.2, 0.3, 0.3), 2)
#' dm.f81 <- dm.setMutationModel(dm, "F84", c(0.3, 0.2, 0.3, 0.2), 2)
#' dm.gtr <- dm.setMutationModel(dm, "GTR", gtr.rates=c(0.2, 0.2, 0.1, 0.1, 0.1, 0.2))
dm.setMutationModel <- function(dm, mutation.model, 
                                base.frequencies, tstv.ratio, 
                                gtr.rates) {

  checkType(mutation.model, c("char", "s"), T, F)
  checkType(base.frequencies, c("num"), F, F)
  checkType(tstv.ratio, c("num", "s"), F, F)
  checkType(gtr.rates, c("num"), F, F)

  if (! mutation.model %in% mutation.models) 
    stop("Allowed values: ", paste(mutation.models, collapse=" "))
  
  mutation.model.nr <- seq(along = mutation.models)[mutation.models == mutation.model]
  dm <- addFeature(dm, "mutation.model", "mutation.model", 
                            fixed.value=mutation.model.nr)

  if ( !missing(tstv.ratio) ) {
    if (!mutation.model %in% c("HKY", "F84"))
      stop("This mutation model does not support a ts/tv ratio")
    dm <- addFeature(dm, "tstv.ratio", "tstv.ratio", fixed.value=tstv.ratio)
  }

  if ( !missing(base.frequencies) ) {
    if ( length(base.frequencies) != 4 ) 
        stop("You must enter frequencies for all 4 bases")
    if (!mutation.model %in% c("HKY", "F84")) 
      stop("This mutation model does not support base frequencies")

    dm <- addFeature(dm, "base.freq.A", "base.freq.A", fixed.value=base.frequencies[1])
    dm <- addFeature(dm, "base.freq.C", "base.freq.C", fixed.value=base.frequencies[2])
    dm <- addFeature(dm, "base.freq.G", "base.freq.G", fixed.value=base.frequencies[3])
    dm <- addFeature(dm, "base.freq.T", "base.freq.T", fixed.value=base.frequencies[4])
  }

  if ( !missing(gtr.rates) ) {
    if ( length(gtr.rates) != 6 ) 
        stop("You must enter rates for all 6 posible substitutions")
    if (!mutation.model %in% c("GTR")) 
      stop("You can specify gtr.rates only with the GTR model")

    dm <- addFeature(dm, "gtr.rate.1", "gtr.rate.1", fixed.value=gtr.rates[1])
    dm <- addFeature(dm, "gtr.rate.2", "gtr.rate.2", fixed.value=gtr.rates[2])
    dm <- addFeature(dm, "gtr.rate.3", "gtr.rate.3", fixed.value=gtr.rates[3])
    dm <- addFeature(dm, "gtr.rate.4", "gtr.rate.4", fixed.value=gtr.rates[4])
    dm <- addFeature(dm, "gtr.rate.5", "gtr.rate.5", fixed.value=gtr.rates[5])
    dm <- addFeature(dm, "gtr.rate.6", "gtr.rate.6", fixed.value=gtr.rates[6])
  }

  return(dm)
}


#-------------------------------------------------------------------
# dm.addMutationRateHeterogenity
#-------------------------------------------------------------------
#' Allows the mutation rate on different sites to vary according to 
#' a Gamma Distribution.
#'
#' This function adds a Gamma distributed rate heterogeneity as implemented 
#' in 'seq-gen' to the model.
#'
#' "The [...] model of rate heterogeneity assigns different rates to different
#' sites according to a gamma distribution (Yang, 1993). The distribution is scaled
#' such that the mean rate for all the sites is 1 but 
#' the user must supply a parameter which describes its shape. A low value for this
#' parameter (<1.0) simulates a large degree of site-specific rate heterogeneity
#' and as this value increases the simulated data becomes more rate-homogeneous.
#' This can be performed as a continuous model, i.e. every site has a different
#' rate sampled from the gamma distribution of the given shape, or as a discrete
#' model, i.e. each site falls into one of N rate categories approximating the
#' gamma distribution. For a review of site-specific rate heterogeneity and its
#' implications for phylogenetic analyses, see Yang (1996)."
#' [From the seq-gen homepage http://bioweb2.pasteur.fr/docs/seq-gen ]
#'
#' The Parameter in this text will be referred to as 'alpha'. Simulation a model
#' with rate heterogeneity requires that 'seq-gen' is installed on your system. 
#' 
#' @param dm  The demographic model to which the rate heterogeneity should be added.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'min.alpha' and 'max.alpha' or 
#'            'fixed.alpha'. It will be named 'new.par.name'
#'            If 'FALSE' the argument 'parameter'
#'            will be evaluated instead.
#' @param min.alpha  If you want to estimate the rate heterogeneity, this will be 
#'            used as the smallest possible value.
#' @param max.alpha  Same as min.growth.rate, but the largest possible value.
#' @param fixed.alpha  If specified, the growth rate will not be estimated,
#'            but assumed to have the given value.
#' @param new.par.name Name for the new parameter.
#'            occurs. See \link{dm.addSpeciationEvent} for more information.
#' @param categories.number If this is set, a fixed number of categories will be
#'            used to model the gamma distribution instead of drawing every parameter
#'            seperately (see text).
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "alpha" will use
#'            a parameter with name alpha that you have previously 
#'            created. You can also use R expression here, i.e. "2*alpha"
#'            or "5*M+2*alpha" (if M is another parameter) will also
#'            work (also the latter does not make much sense).
#' @return    The demographic model with a size change.
#' @export
#' @examples
#' # A model with one smaller population
#' dm <- dm.createThetaTauModel(c(20,37), 88)
#' dm <- dm.setMutationModel(dm, "HKY")
#' dm <- dm.addMutationRateHeterogenity(dm, 0.1, 5, new.par.name="alpha")
dm.addMutationRateHeterogenity <- 
  function(dm, min.alpha, max.alpha, fixed.alpha, par.new=T, 
           new.par.name="alpha", parameter, categories.number) {

  if (par.new) parameter <- new.par.name

  dm <- addFeature(dm, "gamma.rate", parameter, min.alpha, max.alpha,
                   fixed.alpha, par.new, NA, NA, NA)

  if (!missing(categories.number)) {
    dm <- addFeature(dm, "gamma.categories", 
                     fixed.value=categories.number,
                     parameter="gamma.categories")
  }

  return(dm)
}

#-------------------------------------------------------------------
# dm.createThetaTauModel
#-------------------------------------------------------------------
#' Creates a standard "Theta/Tau" demopraphic model.
#' 
#' @param sample.sizes   A numeric vector with the sample sizes of the two pop.sources.
#' @param loci.num       The number of Loci to simulate.
#' @param seq.length     For recombination, each locus is assumed to be of this
#'                       length
#' @return               A Theta/Tau Model
#' @export
#'
#' @examples
#' dm <- dm.createThetaTauModel(c(20,25), 100)
#' dm
dm.createThetaTauModel <- function(sample.sizes, loci.num, seq.length=1000) {
  dm <- dm.createDemographicModel(sample.sizes, loci.num, seq.length)
  dm <- dm.addSpeciationEvent(dm, 0.01, 5, new.time.point.name="tau")
  dm <- dm.addRecombination(dm, fixed.value=20)
  dm <- dm.addMutation(dm, 1, 20)
  return(dm)
}


#-------------------------------------------------------------------
# dm.addOutgroup
#-------------------------------------------------------------------
#' Adds an outgroup to a demographic model
#'
#' This function adds an outgroup consisting of one individual to a 
#' demographic model. The outgroup consists of one individual.
#' An outgroup is required for a finite sites analysis.
#'
#' @param dm The demographic model to which we add the outgroup
#' @param separation.time The time point at which the outgroup splited 
#'           from the ancestral population. This can be an absolute value
#'           (e.g. 10) or relative to another time points (e.g. '5*t_split_1').
#' 
#' @return  The extended demographic model
#' @export
dm.addOutgroup <- function(dm, separation.time) {
  number.of.individuals <- 1
  dm@sampleSizes <- c(dm@sampleSizes, number.of.individuals)
  dm.addSpeciationEvent(dm, in.population=1, 
                        new.time.point=F, 
                        time.point=separation.time) 
}




#-------------------------------------------------------------------
# dm.simSumStats
#-------------------------------------------------------------------
#' Simulates data according to a demographic model
#' 
#' @param dm The demographic model according to which the simulations should be done
#' @param parameters A vector of parameters which should be used for the simulations. 
#'           If a matrix is given, a simulation for each row of the matrix will be performed
#' @param sum.stats A vector with names of the summary statistics to simulate,
#'           or "all" for simulating all possible statistics.
#' @return A matrix where each row is the vector of summary statistics for 
#'         the parameters in the same row of the "parameter" matrix
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm.simSumStats(dm,c(1,10))
dm.simSumStats <- function(dm, parameters, sum.stats=c("all")) {
  checkType(dm, "dm")

  checkParInRange(dm, parameters)

  dm@currentSimProg@singleSimFunc(dm, parameters)
}


scaleDemographicModel <- function(dm, scaling.factor) {
  dm@nLoci <- round(dm@nLoci / scaling.factor)
  return(dm)
}
