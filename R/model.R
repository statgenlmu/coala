#---------------------------------------------------------------
# DemographicModel.R
# Class for representing a model of the evolutionary development
# of two different species.
# 
# Authors:  Paul R. Staab & Lisha Mathew 
# Email:    staab ( at ) bio.lmu.de
# Licence:  GPLv3 or later
#--------------------------------------------------------------

#' @include dm_sim_program.R

setClass("DemographicModel" ,
         representation(features="data.frame",
                        parameters="data.frame",
                        sum.stats="data.frame",
                        tsTvRatio="numeric",
                        finiteSites="logical",
                        currentSimProg="character",
                        options="list",
                        finalized='logical',
                        loci="data.frame")
         )


#-----------------------------------------------------------------------
# Initialization
#-----------------------------------------------------------------------

createFeatureTable <- function(type=character(), parameter=character(),
                               pop.source=numeric(), pop.sink=numeric(),
                               time.point=character(), group=numeric()) {
  
  stopifnot(is.character(type))
  stopifnot(is.character(parameter))
  stopifnot(is.na(pop.source) | is.numeric(pop.source))
  stopifnot(is.na(pop.sink) | is.numeric(pop.sink))
  stopifnot(is.na(time.point) | is.character(time.point))
  stopifnot(is.numeric(group))
  
  data.frame(type=type,
             parameter=parameter,
             pop.source=pop.source,
             pop.sink=pop.sink,
             time.point=time.point,
             group=group,
             stringsAsFactors=F)
}

.init <- function(.Object, sample.size, loci.number, loci.length,
                  finiteSites, tsTvRatio){
  
  .Object <- resetSumStats(.Object)
  
  .Object@features <- createFeatureTable()
  
  .Object@loci <- data.frame(group=numeric(),
                             number=numeric(),  
                             name=character(), 
                             name_l=character(),
                             name_r=character(),
                             length_l=numeric(),
                             length_il=numeric(),
                             length_m=numeric(),
                             length_ir=numeric(),
                             length_r=numeric(),
                             stringsAsFactors=F )

  .Object@parameters <- data.frame(parameter=character(),
                                   lower.range=numeric(),
                                   upper.range=numeric(),
                                   stringsAsFactors=F )  

  for (pop in seq(along = sample.size)) {
    .Object <- .Object + feat_sample(sample.size[pop], pop)
  }
  .Object <- dm.addLocus(.Object, length = loci.length, number = loci.number)
  
  .Object@finiteSites     <- finiteSites
  .Object@tsTvRatio       <- tsTvRatio
  .Object@options         <- list()
  .Object@finalized       <- FALSE

  return(.Object)
}

setMethod("initialize","DemographicModel",.init)
rm(.init)

#-----------------------------------------------------------------------
# Print
#-----------------------------------------------------------------------
.showModel <- function(object) {
  if (!object@finalized) dm = dm.finalize(object)
  .print("Used simulation program:", object@currentSimProg)
  .print()
  
  # Print parameters that get estimated
  
  pars.est = object@parameters
  rownames(pars.est) <- NULL
  
  if (nrow(pars.est) == 0) .print("No Parameters.")
  else { 
    .print("Parameters:")
    print(pars.est)
  }
  .print()
    
  # Print simulation command
  .print("Simulation command:")
  getSimProgram(object@currentSimProg)$print_cmd_func(object)
}

.show <- function(object) {
  object <- dm.finalize(object)
  
  if (is.null(object@options$grp.models)) {
    .showModel(object)
  } else {
    for (group in names(object@options$grp.models)) {
      .print('----------------------------------')
      .print("Group", group)
      .print('----------------------------------')
      .showModel(object@options$grp.models[[group]])
      .print()
    }
  }
} 
setMethod("show", "DemographicModel", .show)
rm(.show)


#------------------------------------------------------------------------------
# Private functions
#------------------------------------------------------------------------------

dm.addParameter <- function(dm, par.name, lower.boundary, upper.boundary) {

  checkType(dm,          c("dm",   "s"), T, F)
  checkType(par.name,    c("char", "s"), T, F)
  checkType(lower.boundary, c("num",  "s"), T, T)
  checkType(upper.boundary, c("num",  "s"), T, T)

  if (par.name %in% dm.getParameters(dm))
    stop("There is already a parameter with name ", par.name)
 
  dm <- appendToParameters(dm, par.name, lower.boundary, upper.boundary)

  dm@finalized = FALSE
  return(dm)
}


dm.addFeature <- function(dm, feature) {
  stopifnot(is.feature(feature))
  
  dm@features <- rbind(dm@features, feature$get_table())
  for (parameter in feature$get_parameters()) {
    dm <- dm + parameter
  }
  
  if (feature$get_inter_locus_var()) {
    dm <- addInterLocusVariation(dm, feature$get_group())
  }
  
  dm
}

#' Adds a summary statistic to the model.
#' 
#' This summary statistic will the be calulated for each simualtion
#' and returned by the simSumStats function.
#' 
#' Avaible summary statistics are
#' 'jsfs' - calculates the Joint Site Frequency Spectrum
#' 'seg.sites' - return the simulated segregating sites as matrix
#' 'file' - returns a file in which the simulation output is written
#' 'fpc' - calculates the Four-Gamete-Condition based statistic
#' 
#' @param dm The demographic model to which a summary statistic should be added. 
#' @param sum.stat The summary statistic to add. Use the names mentioned above.
#' @param group If given, the summary statistic is only calculated for a 
#'        given group of loci.
#' @param population The population for which the summary statistic is calculated.
#'   Currently only used for 'fpc' statistics.
#' @return The model with a summary statistic added.
#' @export
#' @examples
#' dm <- dm.createDemographicModel(c(15, 20), 100)
#' dm <- dm.addSummaryStatistic(dm, 'seg.sites')
dm.addSummaryStatistic <- function(dm, sum.stat, population = 0, group = 0) {
  checkType(dm, "dm")
  checkType(sum.stat, "char")

  if (sum.stat == 'fpc' & population == 0) {
    dm <- dm.addSummaryStatistic(dm, 'fpc', 1, group)
    dm <- dm.addSummaryStatistic(dm, 'fpc', 2, group)
    return(dm)
  }
  
  # Add the summary statistic
  dm@sum.stats = rbind(dm@sum.stats, data.frame(name=sum.stat, 
                                                population = population, 
                                                group=group))
  dm@finalized = FALSE
  
  # Check if there is any simulation program supporting this summary statistic
  for (sim.prog in .jaatha$sim_progs) {
    if (sum.stat %in% sim.prog$possible_sum_stats) return(dm)
  }
  stop("No simulation program for summary statistic", sum.stat)
}



# Helper function that appends a parameter to the "parameters" dataframe
# Does not check the feature for consistency
# This should only be used by addFeature().
appendToParameters <- function(dm, name, lower.range, upper.range) {

  new.parameter <- data.frame(name=name,
                              lower.range=lower.range,
                              upper.range=upper.range,
                              stringsAsFactors=F)
  
  dm@parameters <- rbind(dm@parameters, new.parameter)
  return(dm)
}


# Gets the availible populations
getPopulations <- function(dm){
  unique(searchFeature(dm, 'sample')$pop.source)
}

# Checks if a vector of parameters is within the ranges of the model
checkParInRange <- function(dm, param) {
  if (length(param) != dm.getNPar(dm)) stop("Wrong number of parameters")

  ranges <- dm.getParRanges(dm)
  in.range <- all(ranges[, 1]-1e-11 <= param & param <= ranges[, 2]+1e-11)
  if (!in.range) stop("Parameter combination out of range")
}

# Selects a program for simulation that is capable of all current features
dm.selectSimProg <- function(dm) {
  name <- NULL
  priority <- -Inf
  
  for (sim_prog in .jaatha$sim_progs) {
    if (all(dm@features$type %in% sim_prog$possible_features) & 
        all(dm@sum.stats$name %in% sim_prog$possible_sum_stats)) {
      
      if (sim_prog$priority > priority) {
        name <- sim_prog$name
        priority <- sim_prog$priority
      }
      
    }
  }
  
  if (is.null(name)) stop("No suitable simulation software found!")
  
  dm@currentSimProg <- name
  return(dm)
}

dm.finalize <- function(dm) {
  if (length(dm.getGroups(dm)) == 1) {
    dm <- generateGroupModel(dm, 1)
    dm <- dm.selectSimProg(dm)
    return(getSimProgram(dm@currentSimProg)$finalization_func(dm))
  }

  dm@options$grp.models <- list()
  dm@currentSimProg <- "groups"
  dm.raw <- dm
  
  for (group in dm.getGroups(dm)) {
    grp.model <- generateGroupModel(dm.raw, group)
    grp.model <- dm.finalize(grp.model)
    dm@options$grp.models[[as.character(group)]] <- grp.model
  }
  
  dm@finalized = TRUE
  dm
}



#------------------------------------------------------------------------------
# Getters & Setters for Jaatha
#------------------------------------------------------------------------------
dm.getParameters <- function(dm) {
  checkType(dm,"dm")
  dm@parameters$name
}

dm.getNPar <- function(dm){
  checkType(dm, "dm")
  return(length(dm.getParameters(dm)))
}

dm.getParRanges <- function(dm){
  #checkType(dm,"dm")
  par.ranges <- dm@parameters[ , c("lower.range","upper.range")]
  rownames(par.ranges) <- dm@parameters[ ,"name"]
  par.ranges
}

getThetaName <- function(dm, outer=FALSE, group=0) {
  if (outer) {
    feat <- searchFeature(dm, "mutation_outer", group=group)
    if (nrow(feat) == 0) {
      feat <- searchFeature(dm, "mutation", group=group)
    }
  }  else {
    feat <- searchFeature(dm, "mutation", group=group)
  }
  if (nrow(feat) != 1) stop("Failed to determine mutation rate")
  feat[1, 'parameter']
}

resetSumStats <- function(dm) {
  dm@sum.stats <- dm@sum.stats[FALSE, ]
  dm
}



#------------------------------------------------------------------------------
# Creation new models
#------------------------------------------------------------------------------

#' Create a basic demographic model
#' 
#' This function creates a basic empty demographic model, which
#' is returned. Features like mutation, pop.source splits and 
#' migration can be added afterwards.
#'
#' @param sample.sizes Number of haploid individuals/chromosomes that are sampled. If your model 
#'            consists of multiple populations, this needs to be a vector
#'            containing the sample sizes from each population.
#' @param loci.num     Number of loci that will be simulated
#' @param seq.length   (Average) number of bases for each locus
#' @return            The demographic model
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sample.sizes=c(25,25), loci.num=100)
#' dm <- dm.addSpeciationEvent(dm, 0.01, 5, 'tau', 1, 2)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.createDemographicModel <- function(sample.sizes, loci.num, seq.length=1000) {
  dm <- new("DemographicModel", sample.sizes, loci.num, seq.length, F, .33)
  dm <- dm.addSummaryStatistic(dm, 'jsfs')
  return(dm)
}


# Low level function for adding a locus
addLocus <- function(dm, group=0, number=1, 
                     name='', name_l='', name_r='',
                     length_l=0, length_il=0, length_m=0, 
                     length_ir=0, length_r=0) {
  if (number > 1 & any(dm@loci[dm@loci$group == group, 'number'] != 1)) {
    stop("You can only have multiple loci in one group if 'number' is 1 for all")
  }
  
  dm@loci <- rbind(dm@loci, data.frame(group=group,
                                       number=number,  
                                       name=name, 
                                       name_l=name_l,
                                       name_r=name_r,
                                       length_l=length_l,
                                       length_il=length_il,
                                       length_m=length_m,
                                       length_ir=length_ir,
                                       length_r=length_r))
  
  invisible(dm)
}

#' Defines how many identical loci belong to a group of loci
#'
#' @param dm The Demographic Model
#' @param number The number of loci to add to the group.
#' @param length The average length of the loci or the acctual of the locus 
#'               if just one is used.
#' @param group The group for which we set the loci number
#' @return The changed Demographic Model
#' @export
#' @examples
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addLocus(dm, number = 200, length = 250, group = 1)
dm.addLocus <- function(dm, length, number = 1, group=0) {
  checkType(dm, 'dm')
  if (!is.numeric(number)) stop("'number' needs to be numeric")
  if (!is.numeric(length)) stop("'length' needs to be numeric")
  if (!is.numeric(group)) stop("'group' needs to be numeric")
  
  addLocus(dm, group=group, number=number, length_m=length)
}

#' Adds a trio of loci to a group
#' @param dm The Demographic Model
#' @param locus_names A vector of 3 strings, giving the names for the loci.
#'   The names are used for identifying the loci later (left, middle and right).
#' @param locus_length An integer vector of length 3, giving the length of each 
#'   of the three loci (left, middle and right).
#' @param distance A vector of two, giving the distance between left and middle,
#'   and middle an right locus, in basepairs. 
#' @param group The group to which to add the trio
#' @return The extended demographic model
#' @export
#' @examples
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addLocusTrio(dm, locus_names = c('Solyc00g00500.2',
#'                                           'Solyc00g00520.1',
#'                                           'Solyc00g00540.1'),
#'                       locus_length=c(1250, 1017, 980), 
#'                       distance=c(257, 814))
dm.addLocusTrio <- function(dm, locus_names=c(left='', middle='', right=''), 
                            locus_length=c(left=1000, middle=1000, right=1000), 
                            distance=c(left_middle=500, middle_right=500),
                            group=1) {
  
  checkType(dm, 'dm')
  if (!is.character(locus_names)) stop("'name' needs to be numeric")
  if (length(locus_names) != 3) stop("'name' needs to be a vector of three names")  
  if (!is.numeric(locus_length)) stop("'locus_length' needs to be numeric")
  if (length(locus_length) != 3) 
    stop("'locus_length' needs to be a vector of three names")    
  if (!is.numeric(group)) stop("'group' needs to be numeric")
  
  if (nrow(searchFeature(dm, 'locus_trios', group = group)) == 0) {
    dm <- addFeature(dm, 'locus_trios', parameter = NA, group = group)
  }
  
  addLocus(dm, group=group, 
           name_l = locus_names[1], 
           name = locus_names[2], 
           name_r = locus_names[3],
           length_l=locus_length[1],
           length_il=distance[1],
           length_m=locus_length[2],
           length_ir=distance[2],
           length_r=locus_length[3])
}

# Legacy function for unit testing
dm.setLociNumber <- function(dm, number, group = 0) {
  if (sum(dm@loci$group == group) != 1) stop('More the one set of loci for this group')
  dm@loci[dm@loci$group == group, 'number'] <- number
  dm
}

# Legacy function for unit testing
dm.setLociLength <- function(dm, length, group = 0) {
  if (sum(dm@loci$group == group) != 1) stop('More the one set of loci for this group')
  dm@loci[dm@loci$group == group, 'length_m'] <- length
  dm
}



#' Gets how many loci belong to a group of loci
#'
#' @param dm The Demographic Model
#' @param group The group for which we get the number of loci. Defaults to 
#'              the first group.
#' @return The number of loci in the group
#' @export
#' @examples
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addLocus(dm, number = 200, length = 250, group = 1)
#' dm.getLociNumber(dm)
#' dm.getLociNumber(dm, group = 0)
dm.getLociNumber <- function(dm, group=1) {
  number <- dm@loci[dm@loci$group == group, 'number']
  if (length(number) == 0) {
    if (group > 0) number <- dm@loci[dm@loci$group == 0, 'number']
    else stop("Failed to determine loci number")
  }
  as.integer(sum(number))
}

#' Gets how long the loci in a group are
#'
#' @param dm The Demographic Model
#' @param group The group for which we get the length of loci
#' @return The length of the loci in the group
#' @export
#' @examples
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addLocus(dm, number = 200, length = 250, group = 1)
#' dm.getLociLength(dm)
#' dm.getLociLength(dm, group = 0)
#' dm.getLociLength(dm, group = 2)
dm.getLociLength <- function(dm, group=1) {
  length <- dm@loci[dm@loci$group == group, c(6,8,10), drop = FALSE]
  if (nrow(length) == 0) {
    length <- dm@loci[dm@loci$group == 0, c(6,8,10), drop = FALSE]
  }
  if (nrow(length) == 0) stop("Failed to determine loci length")
  as.integer(sum(length))
}

dm.getLociLengthMatrix <- function(dm, group=1) {
  # Select the rows of the group
  rows <- which(dm@loci$group == group)
  if (sum(rows) == 0) rows <- which(dm@loci$group == 0)
  
  # Repeat the row if number > 1
  if (length(rows) == 1) rows <- rep(rows, dm@loci$number[rows])
  
  # Return the matrix
  llm <- dm@loci[rows, 6:10, drop = FALSE]
  row.names(llm) <- NULL
  as.matrix(llm)
}

dm.getSampleSize <- function(dm) {
  feat.samples <- searchFeature(dm, type="sample")
  stopifnot(nrow(feat.samples) > 0)

  sample.size <- rep(0, max(na.omit(dm@features$pop.source)))
  for (row.nr in 1:nrow(feat.samples)) {
    stopifnot(sample.size[feat.samples$pop.source[row.nr]] == 0)
    sample.size[feat.samples$pop.source[row.nr]] <-
      as.integer(feat.samples$parameter[row.nr])
  }

  sample.size
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
#' fraction of N0, the present day size of population one.
#'
#' If you want to add a slow, continuous change over some time,
#' then use the \link{dm.addGrowth} function.
#' 
#' @param dm  The demographic model to which the size change should be added.
#' @param min.size.factor  If you want to estimate the size factor, this will be 
#'            used as the smallest possible value.
#' @param max.size.factor  Same as min.size.factor, but the largest possible value.
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
#' dm <- dm.createDemographicModel(c(20,37), 88)
#' dm <- dm.addSizeChange(dm, 0.1, 1, population=2, at.time="0")
dm.addSizeChange <- function(dm, min.size.factor, max.size.factor,
                             parameter="q",
                             population, at.time="0") {

  checkType(population, c("num",  "s"), T, F)
  checkType(at.time,    c("char", "s"), T, F)

  dm <- addFeature(dm, "size.change", parameter, min.size.factor, max.size.factor,
                   population, NA, at.time)

  dm
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
#' @param min.growth.rate  If you want to estimate the growth rate, this will be 
#'            used as the smallest possible value.
#' @param max.growth.rate  Same as min.growth.rate, but the largest possible value.
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
#' dm <- dm.createDemographicModel(c(20,37), 88)
#' dm <- dm.addGrowth(dm, 0.1, 2, population=2, at.time="0")
dm.addGrowth <- function(dm, min.growth.rate=NA, max.growth.rate=NA, 
                         parameter="alpha", population, at.time="0") {

  checkType(population, c("num",  "s"), T, F)
  checkType(at.time,    c("char", "s"), T, F)

  dm <- addFeature(dm, "growth", parameter, min.growth.rate, max.growth.rate,
                   population, NA, at.time)
  dm
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
#' dm <- dm.createDemographicModel(10:11, 10, 100)
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

  if (! mutation.model %in% sg.mutation.models) 
    stop("Possible mutation models: ", paste(sg.mutation.models, collapse=" "))
  
  dm <- addFeature(dm, "mutation.model", mutation.model)

  if ( !missing(tstv.ratio) ) {
    if (!mutation.model %in% c("HKY", "F84"))
      stop("This mutation model does not support a ts/tv ratio")
    dm <- addFeature(dm, "tstv.ratio", tstv.ratio)
  }

  if ( !missing(base.frequencies) ) {
    if ( length(base.frequencies) != 4 ) 
        stop("You must enter frequencies for all 4 bases")
    if (!mutation.model %in% c("HKY", "F84")) 
      stop("This mutation model does not support base frequencies")

    dm <- addFeature(dm, "base.freq.A", base.frequencies[1])
    dm <- addFeature(dm, "base.freq.C", base.frequencies[2])
    dm <- addFeature(dm, "base.freq.G", base.frequencies[3])
    dm <- addFeature(dm, "base.freq.T", base.frequencies[4])
  }

  if ( !missing(gtr.rates) ) {
    if ( length(gtr.rates) != 6 ) 
        stop("You must enter rates for all 6 posible substitutions")
    if (!mutation.model %in% c("GTR")) 
      stop("You can specify gtr.rates only with the GTR model")

    dm <- addFeature(dm, "gtr.rate.1", gtr.rates[1])
    dm <- addFeature(dm, "gtr.rate.2", gtr.rates[2])
    dm <- addFeature(dm, "gtr.rate.3", gtr.rates[3])
    dm <- addFeature(dm, "gtr.rate.4", gtr.rates[4])
    dm <- addFeature(dm, "gtr.rate.5", gtr.rates[5])
    dm <- addFeature(dm, "gtr.rate.6", gtr.rates[6])
  }

  return(dm)
}


#-------------------------------------------------------------------
# dm.addMutationRateHeterogenity
#-------------------------------------------------------------------
#' Allows the mutation rate on different sites within one locus to 
#' vary according to a Gamma Distribution.
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
#' @param min.alpha  If you want to estimate the rate heterogeneity, this will be 
#'            used as the smallest possible value.
#' @param max.alpha  Same as min.growth.rate, but the largest possible value.
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
#' dm <- dm.createDemographicModel(c(20,37), 88)
#' dm <- dm.setMutationModel(dm, "HKY")
#' dm <- dm.addMutationRateHeterogenity(dm, 0.1, 5, parameter="alpha")
dm.addMutationRateHeterogenity <- 
  function(dm, min.alpha, max.alpha, parameter="alpha", categories.number) {

  dm <- addFeature(dm, "gamma.rate", parameter, min.alpha, max.alpha, NA, NA, NA)

  if (!missing(categories.number)) {
    dm <- addFeature(dm, "gamma.categories", parameter=categories.number)
  }

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
#' @param separation_time The time point at which the outgroup splits 
#'           from the ancestral population. This can be an absolute value
#'           (e.g. 10) or relative to another time points (e.g. '5*t_split_1').
#' @param sample_size The number of individuals in the group
#' @param anc_pop The ancestral population of the non-outgroup individuals. The
#'          outgroup splits from this population at the separation_time.
#' @return  The extended demographic model
#' @export
dm.addOutgroup <- function(dm, separation_time, sample_size = 1, anc_pop = 1) {
  pop <- max(na.omit(dm@features$pop.source)) + 1
  dm <- addFeature(dm, "sample", as.character(sample_size), 
                     pop.source=pop, group=0, time.point='0')
  dm <- dm.addSpeciationEvent(dm, in.pop=anc_pop, to.pop=pop,
                              time.point=separation_time)
  dm <- addFeature(dm, "outgroup",sample_size, 
                   pop.source = anc_pop, pop.sink = pop)
  dm
}

dm.getOutgroupSize <- function(dm) {
  as.integer(searchFeature(dm, 'outgroup')$parameter)
} 


#' Adds positiv selection to a model
#' 
#' @inheritParams dm.addMutation
#' @param min.strength Minimal strength of selection
#' @param max.strength Maximal strength of selection
#' @param fraction.neutral Optionally, a fraction of the loci in the group can
#'   be neutral.
#' @param population The populaton in which the allele is selected.
#' @param at.time The time at which the selection starts.
#' @export
dm.addPositiveSelection <- function(dm, min.strength=NA, max.strength=NA,
                         parameter='s', variance = 0, fraction.neutral = 0,
                         population, at.time, group=0) {

  checkType(population, c("num",  "s"), T, F)
  checkType(at.time,    c("char", "s"), T, F)

  dm <- addFeature(dm, "pos.selection", parameter, min.strength, max.strength,
                   population, NA, at.time, group,
                   variance = variance, zero.inflation = fraction.neutral)

  dm
}


#' Adds balancing selection to a model
#' 
#' @inheritParams dm.addMutation
#' @param min.strength Minimal strength of selection
#' @param max.strength Maximal strength of selection
#' @param fraction.neutral Optionally, a fraction of the loci in the group can
#'   be neutral.
#' @param population The populaton in which the allele is selected.
#' @param at.time The time at which the selection starts.
#' @export
dm.addBalancingSelection <- function(dm, min.strength=NA, max.strength=NA,
                                    parameter='s', variance = 0, fraction.neutral = 0,
                                    population, at.time, group=0) {
  
  checkType(population, c("num",  "s"), T, F)
  checkType(at.time,    c("char", "s"), T, F)
  
  dm <- addFeature(dm, "bal.selection", parameter, min.strength, max.strength,
                   population, NA, at.time, group,
                   variance = variance, zero.inflation = fraction.neutral)
  
  dm
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
#' dm <- dm.addSpeciationEvent(dm, 0.01, 5, 'tau', 1, 2)
#' dm <- dm.addMutation(dm,1,20)
#' dm.simSumStats(dm,c(1,10))
dm.simSumStats <- function(dm, parameters, sum.stats=c("all")) {
  checkType(dm, "dm")
  checkParInRange(dm, parameters)
  
  if (!dm@finalized) dm = dm.finalize(dm)

  if (dm@currentSimProg != "groups") {
    return(getSimProgram(dm@currentSimProg)$sim_func(dm, parameters))
  } 
    
  sum.stats <- list(pars=parameters)
  for (group in dm.getGroups(dm)) {
    dm.grp <- dm@options$grp.models[[as.character(group)]]
    sum.stats.grp <- getSimProgram(dm.grp@currentSimProg)$sim_func(dm.grp, parameters)
    for (i in seq(along = sum.stats.grp)) {
      if (names(sum.stats.grp)[i] == 'pars') next()
      name <- paste(names(sum.stats.grp)[i], group, sep='.')
      sum.stats[[name]] <- sum.stats.grp[[i]]
    }
  }
  
  sum.stats
}


generateGroupModel <- function(dm, group) {
  if (all(dm@features$group == 0) & 
      all(dm@sum.stats$group == 0) &
      all(dm@loci$group == 0) ) return(dm)
  
  if (!is.null(dm@options$grp.models[[as.character(group)]])) { 
    return(dm@options$grp.models[[as.character(group)]]) 
  }
  
  # Features
  dm@features <- searchFeature(dm, group = group)
  dm@features$group <- 0
  
  # Sum.Stats
  dm@sum.stats <- dm@sum.stats[dm@sum.stats$group %in% c(0, group), ]
  dm@sum.stats$group <- 0
  
  # Loci
  loci <- dm@loci[dm@loci$group == group, , drop=FALSE]
  if (nrow(loci) > 0) dm@loci <- loci
  else dm@loci <- dm@loci[dm@loci$group == 0, , drop=FALSE]
  dm@loci$group <- 0
  
  # Options
  group.name <- paste("group", group, sep='.')
  if (!is.null(dm@options[[group.name]])) {
    for (option in names(dm@options[[group.name]])) {
      dm@options[[option]] <- dm@options[[group.name]][[option]]
    }
  }
  
  dm
}

searchFeature <- function(dm, type=NULL, parameter=NULL, pop.source=NULL,
                       pop.sink=NULL, time.point=NULL, group=NULL) {

  mask <- rep(TRUE, nrow(dm@features))

  if (!is.null(type)) mask <- mask & dm@features$type %in% type

  if (!is.null(parameter)) {
    if (is.na(parameter)) { 
      mask <- mask & is.na(dm@features$parameter) 
    } else {
      mask <- mask & dm@features$parameter %in% parameter 
    }
  }

  if (!is.null(pop.source)) {
    if (is.na(pop.source)) { 
      mask <- mask & is.na(dm@features$pop.source) 
    } else {
      mask <- mask & dm@features$pop.source %in% pop.source
    }
  }

  if (!is.null(pop.sink)) {
    if (is.na(pop.sink)) { 
      mask <- mask & is.na(dm@features$pop.sink) 
    } else {
      mask <- mask & dm@features$pop.sink %in% pop.sink
    }
  }

  if (!is.null(time.point)) {
    if (is.na(time.point)) { 
      mask <- mask & is.na(dm@features$time.point) 
    } else {
      mask <- mask & dm@features$time.point %in% time.point
    }
  }

  if (!is.null(group)) {
    if (group == 0) mask <- mask & dm@features$group == 0
    else {
      mask <- mask & dm@features$group %in% c(0, group)
      
      # Check if values for the default group are overwritten in the focal group
      grp_0 <- dm@features$group == 0
      overwritten <- grp_0[mask]
      for (i in which(overwritten)) {
        duplicats <- searchFeature(dm, type=dm@features$type[mask][i],
                                   pop.source=dm@features$pop.source[mask][i],
                                   pop.sink=dm@features$pop.sink[mask][i])
        if (sum(duplicats$group == group) == 0) overwritten[i] <- FALSE  
      }
      mask[mask] <- !overwritten
    }
  }
  
  return(dm@features[mask, ])
}

#-------------------------------------------------------------------
# dm.getGroups
#-------------------------------------------------------------------
#' Returns the groups currently in the model
#' 
#' @param dm The demographic model
#' @return The groups in the model.
#' @export
dm.getGroups <- function(dm) {
  if (all(c(dm@features$group == 0, 
            dm@sum.stats$group == 0, 
            dm@loci$group == 0))) return(1)

  groups <- sort(unique(c(1, 
                          dm@features$group, 
                          dm@sum.stats$group, 
                          dm@loci$group)))
  
  return(groups[groups != 0])
}


dm.getSummaryStatistics <- function(dm, group = 1, pop) {
  rows <- dm@sum.stats$group %in% c(0,group)
  if (!missing(pop)) rows <- rows & dm@sum.stats$population == pop
  unique(dm@sum.stats[rows, 'name'])
}


scaleDemographicModel <- function(dm, scaling.factor) {
  for (group in unique(dm@loci$group)) {
    dm <- dm.setLociNumber(dm, 
                     round(dm.getLociNumber(dm, group) / scaling.factor),
                     group)
  }
  dm
}

addInterLocusVariation <- function(dm, group = 0) {
  stopifnot(is.numeric(group))
  if (hasInterLocusVariation(dm, group)) return(dm)
  dm + Feature$new('inter_locus_variation', par_const(NA), group = group)
}

hasInterLocusVariation <- function(dm, group = 0) {
  nrow(searchFeature(dm, 'inter_locus_variation', group = group)) > 0
}

getIndOfPop <- function(dm, population) {
  sasi <- dm.getSampleSize(dm)
  if (population == 1) return(1:sasi[1])
  else if (population == 2) return(1:sasi[2]+sasi[1])
  else stop("Invalid population")
}

#' Set the mutation rates for trios
#' @inheritParams dm.addMutation
#' @param middle_rate The mutation rate used for the middle locus
#' @param outer_rate The mutation rate for the two outer loci
#' @export
dm.setTrioMutationRates <- function(dm, middle_rate, outer_rate, group = 0) {
  dm <- addFeature(dm, 'mutation', parameter = middle_rate, group = group)
  dm <- addFeature(dm, 'mutation_outer', parameter = outer_rate, group = group)
}

dm.hasTrios <- function(dm, group=0) {
  sum(dm.getLociLengthMatrix(dm, group)[,-3]) > 0
}