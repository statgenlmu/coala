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
CoalModel <- function(sample_size=0, loci_number=0, loci_length=1000) {
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

  model$id <- get_id()

  # Add sample sizes
  for (pop in seq(along = sample_size)) {
    if (sample_size[pop] > 0) model <- model + feat_sample(sample_size[pop], pop)
  }

  # Add locus
  if (loci_number > 0) {
    model <- model + locus_averaged(loci_number, loci_length)
  }

  model
}


is.model <- function(model) {
  "CoalModel" %in% class(model)
}


#-----------------------------------------------------------------------
# Print
#-----------------------------------------------------------------------
.showModel <- function(object) {
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
determine_sim_prog <- function(dm) {
  name <- read_cache(dm, 'sim_prog')

  if (is.null(name)) {
    message('Determining simulation program')

    if (length(get_groups(dm)) > 1) {
      name <- 'groups'
    } else {
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
    }

    if (is.null(name)) stop("No suitable simulation software found!")
    cache(dm, 'sim_prog', name)
  }

  name
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


generateGroupModel <- function(dm, group) {
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

  dm$id <- get_id()
  dm
}


get_group_model <- function(model, group) {
  grp_model <- read_cache(model, paste0('grp_model_', group))
  if (is.null(grp_model)) {
    message("Generating group model ", group)
    grp_model <- generateGroupModel(model, group)
    cache(model, paste0('grp_model_', group), grp_model)
  }
  grp_model
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


has_trios <- function(dm, group=0) {
  sum(get_locus_length_matrix(dm, group)[,-3]) > 0
}
