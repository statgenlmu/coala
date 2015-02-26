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
                           name=character(), name_l=character(),
                           name_r=character(),
                           length_l=numeric(), length_il=numeric(),
                           length_m=numeric(),
                           length_ir=numeric(), length_r=numeric(),
                           stringsAsFactors=F )

  model$parameters <- data.frame(parameter=character(),
                                 lower.range=numeric(),
                                 upper.range=numeric(),
                                 stringsAsFactors=F)

  model$sum_stats <- create_sumstat_container()

  model$id <- get_id()

  # Add sample sizes
  for (pop in seq(along = sample_size)) {
    if (sample_size[pop] > 0) model <- model +
      feat_sample(sample_size[pop], pop)
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


# Checks if a vector of parameters is within the ranges of the model
checkParInRange <- function(dm, param) {
  if (length(param) != nrow(get_parameter_table(dm))) {
    stop("Wrong number of parameters")
  }

  ranges <- get_parameter_table(dm)[,2:3]
  in.range <- all(ranges[, 1] - 1e-11 <= param & param <= ranges[, 2] + 1e-11)
  if (!in.range) stop("Parameter combination out of range")
}

# Selects a program for simulation that is capable of all current features
determine_sim_prog <- function(dm) {
  name <- read_cache(dm, 'sim_prog')

  if (is.null(name)) {
    if (length(get_groups(dm)) > 1) {
      name <- 'groups'
    } else {
      priority <- -Inf

      for (sim_prog_name in ls(sim_programs)) {
        sim_prog <- getSimProgram(sim_prog_name)
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


generateGroupModel <- function(dm, group) {
  # Features
  dm$features <- searchFeature(dm, group = group)
  dm$features$group <- 0

  # sum_stats
  dm$sum_stats <- get_group_statistics(dm, group)

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
    grp_model <- generateGroupModel(model, group)
    cache(model, paste0('grp_model_', group), grp_model)
  }
  grp_model
}


searchFeature <- function(dm, type=NULL, pop.source=NULL,
                          pop.sink=NULL, time.point=NULL, group=NULL) {

  mask <- rep(TRUE, nrow(get_feature_table(dm)))

  if (!is.null(type)) mask <- mask & dm$features$type %in% type

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
  if (has_inter_locus_var(dm, group)) return(dm)
  dm + Feature$new('inter_locus_variation', par_const(NA), group = group)
}

has_inter_locus_var <- function(dm, group = 0) {
  nrow(searchFeature(dm, 'inter_locus_variation', group = group)) > 0
}


has_trios <- function(dm, group=0) {
  sum(get_locus_length_matrix(dm, group)[,-3]) > 0
}


# Converts a position on the middle locus to the relative position
# on the simulated stretch
conv_middle_to_trio_pos <- function(pos, model, group=1,
                                    relative_out=TRUE, relative_in=TRUE) {
  llm <- get_locus_length_matrix(model, group)

  pos <- ifelse(relative_in, pos * llm[,3], pos) + llm[,1] + llm[,2]
  if (relative_out) pos <- pos / rowSums(llm)

  pos
}


get_snp_positions <- function(seg_sites, model, group=1, relative=TRUE) {
  assert_that(length(seg_sites) == get_locus_number(model, group))
  llm <- get_locus_length_matrix(model, group)
  lapply(1:length(seg_sites), function(locus) {
    pos <- attr(seg_sites[[locus]], 'position')
    trio_locus <- attr(seg_sites[[locus]], 'locus')
    if (is.null(trio_locus)) trio_locus <- 0
    pos[trio_locus == -1] <- pos[trio_locus == -1] * llm[locus, 1]
    pos[trio_locus == 0] <- pos[trio_locus == 0] * llm[locus, 3] +
      sum(llm[locus, 1:2])
    pos[trio_locus == 1] <- pos[trio_locus == 1] * llm[locus, 5] +
      sum(llm[locus, 1:4])
    if (relative) pos <- pos / sum(llm[locus,])
    pos
  })
}
