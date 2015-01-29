#---------------------------------------------------------------
# DemographicModel.R
# Class for representing a model of the evolutionary development
# of two different species.
#
# Authors:  Paul R. Staab & Lisha Mathew
# Email:    staab ( at ) bio.lmu.de
# Licence:  GPLv3 or later
#--------------------------------------------------------------

#' @include sim_program.R

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

#' @export
CoalModel <- function(sample_size, loci_number, loci_length=1000) {
  model <- list()
  class(model) <- c("CoalModel", class(Base_Object))

  model$features <- createFeatureTable()

  model$loci <- data.frame(group=numeric(), number=numeric(),
                           name=character(), name_l=character(), name_r=character(),
                           length_l=numeric(), length_il=numeric(),
                           length_m=numeric(),
                           length_ir=numeric(), length_r=numeric(),
                           stringsAsFactors=F )

  model$parameters <- data.frame(parameter=character(),
                                 lower.range=numeric(),
                                 upper.range=numeric(),
                                 stringsAsFactors=F)

  model$sum_stats <- data.frame(name=character(),
                                population=numeric(),
                                group=numeric(),
                                stringsAsFactors=F)

  # Add sample sizes
  for (pop in seq(along = sample_size)) {
    model <- model + feat_sample(sample_size[pop], pop)
  }

  # Add loci
  model <- dm.addLocus(model, length = loci_length, number = loci_number)

  model$options <- list()
  model$finalized <- FALSE

  model
}


is.model <- function(model) {
  "CoalModel" %in% class(model)
}


#-----------------------------------------------------------------------
# Print
#-----------------------------------------------------------------------
.showModel <- function(object) {
  if (!object@finalized) dm = dm.finalize(object)
  cat("Used simulation program:", object$currentSimProg, "\n\n")

  # Print parameters that get estimated

  pars.est = object@parameters
  rownames(pars.est) <- NULL

  if (nrow(pars.est) == 0) cat("No Parameters.\n")
  else {
    cat("Parameters:\n")
    print(pars.est)
  }
  cat('\n')

  # Print simulation command
  cat("Simulation command:\n")
  getSimProgram(object$currentSimProg)$print_cmd_func(object)
}


.show <- function(object) {
  object <- dm.finalize(object)

  if (is.null(object@options$grp.models)) {
    .showModel(object)
  } else {
    for (group in names(object@options$grp.models)) {
      cat('----------------------------------\n')
      cat("Group", group, '\n')
      cat('----------------------------------\n')
      .showModel(object@options$grp.models[[group]])
      cat('\n')
    }
  }
}


# Checks if a vector of parameters is within the ranges of the model
checkParInRange <- function(dm, param) {
  if (length(param) != nrow(get_parameter_table(dm))) stop("Wrong number of parameters")

  ranges <- get_parameter_table(dm)[,2:3]
  in.range <- all(ranges[, 1]-1e-11 <= param & param <= ranges[, 2]+1e-11)
  if (!in.range) stop("Parameter combination out of range")
}

# Selects a program for simulation that is capable of all current features
dm.selectSimProg <- function(dm) {
  name <- NULL
  priority <- -Inf

  for (sim_prog_name in ls(sim_programs)) {
    sim_prog = getSimProgram(sim_prog_name)
    if (all(dm$features$type %in% sim_prog$possible_features) &
        all(dm$sum_stats$name %in% sim_prog$possible_sum_stats)) {

      if (sim_prog$priority > priority) {
        name <- sim_prog$name
        priority <- sim_prog$priority
      }

    }
  }

  if (is.null(name)) stop("No suitable simulation software found!")

  dm$currentSimProg <- name
  return(dm)
}

dm.finalize <- function(dm) {
  if (length(get_summary_statistics(dm)) == 0) {
    stop("Model has no summary statistics!")
  }
  if (length(get_groups(dm)) == 1) {
    dm <- generateGroupModel(dm, 1)
    dm <- dm.selectSimProg(dm)
    return(getSimProgram(dm$currentSimProg)$finalization_func(dm))
  }

  dm$options$grp.models <- list()
  dm$currentSimProg <- "groups"
  dm.raw <- dm

  for (group in get_groups(dm)) {
    grp.model <- generateGroupModel(dm.raw, group)
    grp.model <- dm.finalize(grp.model)
    dm$options$grp.models[[as.character(group)]] <- grp.model
  }

  dm$finalized = TRUE
  dm
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
  dm$sum_stats <- dm$sum_stats[FALSE, ]
  dm
}


# Low level function for adding a locus
addLocus <- function(dm, group=0, number=1,
                     name='', name_l='', name_r='',
                     length_l=0, length_il=0, length_m=0,
                     length_ir=0, length_r=0) {

  if (number > 1 & any(dm$loci[dm$loci$group == group, 'number'] != 1)) {
    stop("You can only have multiple loci in one group if 'number' is 1 for all")
  }

  dm$loci <- rbind(dm$loci, data.frame(group=group,
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
#' dm <- CoalModel(c(25,25), 100)
#' dm <- dm.addLocus(dm, number = 200, length = 250, group = 1)
dm.addLocus <- function(dm, length, number = 1, group=0) {
  stopifnot(is.model(dm))
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
#' dm <- CoalModel(c(25,25), 100)
#' dm <- dm.addLocusTrio(dm, locus_names = c('Solyc00g00500.2',
#'                                           'Solyc00g00520.1',
#'                                           'Solyc00g00540.1'),
#'                       locus_length=c(1250, 1017, 980),
#'                       distance=c(257, 814))
dm.addLocusTrio <- function(dm, locus_names=c(left='', middle='', right=''),
                            locus_length=c(left=1000, middle=1000, right=1000),
                            distance=c(left_middle=500, middle_right=500),
                            group=1) {

  stopifnot(is.model(dm))
  if (!is.character(locus_names)) stop("'name' needs to be numeric")
  if (length(locus_names) != 3) stop("'name' needs to be a vector of three names")
  if (!is.numeric(locus_length)) stop("'locus_length' needs to be numeric")
  if (length(locus_length) != 3)
    stop("'locus_length' needs to be a vector of three names")
  if (!is.numeric(group)) stop("'group' needs to be numeric")

  if (nrow(searchFeature(dm, 'locus_trios', group = group)) == 0) {
    dm <- dm + Feature$new('locus_trios', par_const(NA), group = group)
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
  if (sum(dm$loci$group == group) != 1) stop('More the one set of loci for this group')
  dm$loci[dm$loci$group == group, 'number'] <- number
  dm
}

# Legacy function for unit testing
dm.setLociLength <- function(dm, length, group = 0) {
  if (sum(dm$loci$group == group) != 1) stop('More the one set of loci for this group')
  dm$loci[dm$loci$group == group, 'length_m'] <- length
  dm
}


generateGroupModel <- function(dm, group) {
  if (all(dm$features$group == 0) &
      all(dm$sum_stats$group == 0) &
      all(dm$loci$group == 0) ) return(dm)

  if (!is.null(dm$options$grp.models[[as.character(group)]])) {
    return(dm$options$grp.models[[as.character(group)]])
  }

  # Features
  dm$features <- searchFeature(dm, group = group)
  dm$features$group <- 0

  # sum_stats
  dm$sum_stats <- dm$sum_stats[dm$sum_stats$group %in% c(0, group), , drop=FALSE]
  dm$sum_stats$group <- 0

  # Loci
  loci <- dm$loci[dm$loci$group == group, , drop=FALSE]
  if (nrow(loci) > 0) dm$loci <- loci
  else dm$loci <- dm$loci[dm$loci$group == 0, , drop=FALSE]
  dm$loci$group <- 0

  # Options
  group.name <- paste("group", group, sep='.')
  if (!is.null(dm$options[[group.name]])) {
    for (option in names(dm$options[[group.name]])) {
      dm$options[[option]] <- dm$options[[group.name]][[option]]
    }
  }

  dm
}

searchFeature <- function(dm, type=NULL, parameter=NULL, pop.source=NULL,
                          pop.sink=NULL, time.point=NULL, group=NULL) {

  mask <- rep(TRUE, nrow(get_feature_table(dm)))

  if (!is.null(type)) mask <- mask & dm$features$type %in% type

  if (!is.null(parameter)) {
    if (is.na(parameter)) {
      mask <- mask & is.na(dm$features$parameter)
    } else {
      mask <- mask & dm$features$parameter %in% parameter
    }
  }

  if (!is.null(pop.source)) {
    if (is.na(pop.source)) {
      mask <- mask & is.na(dm$features$pop.source)
    } else {
      mask <- mask & dm$features$pop.source %in% pop.source
    }
  }

  if (!is.null(pop.sink)) {
    if (is.na(pop.sink)) {
      mask <- mask & is.na(dm$features$pop.sink)
    } else {
      mask <- mask & dm$features$pop.sink %in% pop.sink
    }
  }

  if (!is.null(time.point)) {
    if (is.na(time.point)) {
      mask <- mask & is.na(dm$features$time.point)
    } else {
      mask <- mask & dm$features$time.point %in% time.point
    }
  }

  if (!is.null(group)) {
    if (group == 0) mask <- mask & dm$features$group == 0
    else {
      mask <- mask & dm$features$group %in% c(0, group)

      # Check if values for the default group are overwritten in the focal group
      grp_0 <- dm$features$group == 0
      overwritten <- grp_0[mask]
      for (i in which(overwritten)) {
        duplicats <- searchFeature(dm, type=dm$features$type[mask][i],
                                   pop.source=dm$features$pop.source[mask][i],
                                   pop.sink=dm$features$pop.sink[mask][i])
        if (sum(duplicats$group == group) == 0) overwritten[i] <- FALSE
      }
      mask[mask] <- !overwritten
    }
  }

  dm$features[mask, ]
}




scaleDemographicModel <- function(dm, scaling.factor) {
  for (group in unique(dm$loci$group)) {
    dm <- dm.setLociNumber(dm,
                     round(get_locus_number(dm, group) / scaling.factor),
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

#' Set the mutation rates for trios
#' @param middle_rate The mutation rate used for the middle locus
#' @param outer_rate The mutation rate for the two outer loci
#' @export
dm.setTrioMutationRates <- function(dm, middle_rate, outer_rate, group = 0) {
  dm <- addFeature(dm, 'mutation', parameter = middle_rate, group = group)
  dm <- addFeature(dm, 'mutation_outer', parameter = outer_rate, group = group)
}

dm.hasTrios <- function(dm, group=0) {
  sum(get_locus_length_matrix(dm, group)[,-3]) > 0
}
