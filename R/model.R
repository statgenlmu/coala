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

create_feature_table <- function(type=character(), parameter=character(),
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
coal_model <- function(sample_size=0, loci_number=0, loci_length=1000) {
  model <- list()
  class(model) <- c("coal_model")

  model$feature_table <- create_feature_table()
  model$features <- list()

  model$loci <- data.frame(group=numeric(), number=numeric(),
                           name=character(), name_l=character(),
                           name_r=character(),
                           length_l=numeric(), length_il=numeric(),
                           length_m=numeric(),
                           length_ir=numeric(), length_r=numeric(),
                           stringsAsFactors=F )

  model$parameter <- list()
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
  "coal_model" %in% class(model)
}


# Checks if a vector of parameters is within the ranges of the model
check_par_range <- function(dm, param) {
  if (length(param) != nrow(get_parameter_table(dm))) {
    stop("Wrong number of parameters")
  }

  ranges <- get_parameter_table(dm)[,2:3]
  in.range <- all(ranges[, 1] - 1e-11 <= param & param <= ranges[, 2] + 1e-11)
  if (!in.range) stop("Parameter combination out of range")
}

# Selects a program for simulation that is capable of all current features
determine_simprog <- function(dm) {
  name <- read_cache(dm, 'simprog')

  if (is.null(name)) {
    if (length(get_groups(dm)) > 1) {
      name <- 'groups'
    } else {
      priority <- -Inf

      for (simprog_name in ls(simprograms)) {
        simprog <- get_simprog(simprog_name)
        if (all(get_feature_table(dm)$type %in% simprog$possible_features) &
              all(dm$sum_stats$name %in% simprog$possible_sum_stats)) {

          if (simprog$priority > priority) {
            name <- simprog$name
            priority <- simprog$priority
          }
        }
      }
    }

    if (is.null(name)) stop("No suitable simulation software found!")
    cache(dm, 'simprog', name)
  }

  name
}


get_mutation_par <- function(dm, outer=FALSE, group=0) {
  if (outer) {
    feat <- search_feature(dm, "mutation_outer", group=group)
    if (nrow(feat) == 0) {
      feat <- search_feature(dm, "mutation", group=group)
    }
  }  else {
    feat <- search_feature(dm, "mutation", group=group)
  }
  if (nrow(feat) != 1) stop("Failed to determine mutation rate")
  feat[1, 'parameter']
}


create_group_model <- function(dm, group) {
  # Features
  dm$feature_table <- search_feature(dm, group = group)
  dm$feature_table$group <- 0

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
    grp_model <- create_group_model(model, group)
    cache(model, paste0('grp_model_', group), grp_model)
  }
  grp_model
}


search_feature <- function(dm, type=NULL, pop.source=NULL,
                          pop.sink=NULL, time.point=NULL, group=NULL) {

  feat_tbl <- get_feature_table(dm)
  mask <- rep(TRUE, nrow(feat_tbl))

  if (!is.null(type)) mask <- mask & feat_tbl$type %in% type

  if (!is.null(pop.source)) {
    if (is.na(pop.source)) {
      mask <- mask & is.na(feat_tbl$pop.source)
    } else {
      mask <- mask & feat_tbl$pop.source %in% pop.source
    }
  }

  if (!is.null(pop.sink)) {
    if (is.na(pop.sink)) {
      mask <- mask & is.na(feat_tbl$pop.sink)
    } else {
      mask <- mask & feat_tbl$pop.sink %in% pop.sink
    }
  }

  if (!is.null(time.point)) {
    if (is.na(time.point)) {
      mask <- mask & is.na(feat_tbl$time.point)
    } else {
      mask <- mask & feat_tbl$time.point %in% time.point
    }
  }

  if (!is.null(group)) {
    if (group == 0) mask <- mask & feat_tbl$group == 0
    else {
      mask <- mask & feat_tbl$group %in% c(0, group)

      # Check if values for the default group are overwritten in the focal group
      grp_0 <- feat_tbl$group == 0
      overwritten <- grp_0[mask]
      for (i in which(overwritten)) {
        duplicats <- search_feature(dm, type=feat_tbl$type[mask][i],
                                    pop.source=feat_tbl$pop.source[mask][i],
                                    pop.sink=feat_tbl$pop.sink[mask][i])
        if (sum(duplicats$group == group) == 0) overwritten[i] <- FALSE
      }
      mask[mask] <- !overwritten
    }
  }

  feat_tbl[mask, ]
}


add_inter_locus_var <- function(dm, group = 0) {
  stopifnot(is.numeric(group))
  if (has_inter_locus_var(dm, group)) return(dm)
  dm + Feature$new('inter_locus_variation', par_const(NA), group = group)
}

has_inter_locus_var <- function(dm, group = 0) {
  nrow(search_feature(dm, 'inter_locus_variation', group = group)) > 0
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
