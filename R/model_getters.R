#' Getters for coalescent models
#'
#' @param model The coalescent model from which aspects are returned
#'
#' @export
#' @author Paul Staab
#' @describeIn get_feature_table Returns the features of a model as a data.frame
get_feature_table <- function(model) {
  stopifnot(is.model(model))

  feat_table <- read_cache(model, "feature_table")
  if (is.null(feat_table)) {
    if (length(get_features(model)) == 0) {
      feat_table <- create_feature_table()
    } else {
      feat_table <- do.call(rbind, lapply(get_features(model), function (feat) {
        feat$get_table()
      }))
    }
    cache(model, "feature_table", feat_table)
  }

  feat_table
}


get_features <- function(model) model$features


#' @export
#' @describeIn get_feature_table Returns the ranged parameters of a model as a
#'   data.frame
get_parameter_table <- function(model) {
  stopifnot(is.model(model))

  par_table <- read_cache(model, "par_table")
  if (is.null(par_table)) {
    if (length(get_parameter(model)) == 0) {
      par_table <- (data.frame(name=character(),
                               lower.range=numeric(),
                               upper.range=numeric(),
                               stringsAsFactors=F))
    } else {
      par_table <- do.call(rbind, lapply(get_parameter(model), function (par) {
        data.frame(name=par$get_name(),
                   lower.range=par$get_range()[1],
                   upper.range=par$get_range()[2],
                   stringsAsFactors=F)
      }))
    }
    cache(model, "par_table", par_table)
  }

  par_table
}


#' @export
#' @describeIn get_feature_table Returns the ranged parameters of a model
get_parameter <- function(model) {
  stopifnot(is.model(model))
  model$parameter
}


#' @param group The locus group for which aspects are returned
#' @describeIn get_feature_table Returns the length of the loci in a locus group
get_locus_length <- function(dm) {
  ll <- read_cache(dm, "locus_length")

  if (is.null(ll)) {
    ll <- as.integer(rowSums(
      get_locus_length_matrix(dm, FALSE)[, c(1, 3, 5), drop = FALSE]
    ))

    cache(dm, "locus_length", ll)
  }

  ll
}


#' @describeIn get_feature_table Returns the number of loci in a locus group
#' @export
get_locus_number <- function(dm) {
  number <- sum(sapply(dm$loci, function(x) {
    n <- x$get_number()
    if (n > 1) n <- n / dm$scaling_factor
    n
  }))

  max(as.integer(round(number)))
}


#' @describeIn get_feature_table Returns a vector of populations in the model
#' @export
get_populations <- function(dm) {
  unique(search_feature(dm, 'sample')$pop.source)
}


#' @describeIn get_feature_table Returns a vector of samples sizes per
#'   population.
#' @param for_sim If true, the sample size used internally for the simulation
#'   will be reported rather than the number of actuall samples. The numbers
#'   can be unequal for the simulation of unphased data.
#' @export
get_sample_size <- function(dm, for_sim=FALSE) {
  feat.samples <- search_feature(dm, type="sample")
  stopifnot(nrow(feat.samples) > 0)

  sample_size <- rep(0, length(get_populations(dm)))
  for (row.nr in 1:nrow(feat.samples)) {
    stopifnot(sample_size[feat.samples$pop.source[row.nr]] == 0)
    sample_size[feat.samples$pop.source[row.nr]] <-
      as.integer(feat.samples$parameter[row.nr])
  }

  if (for_sim) {
    sample_size <- sample_size * get_ploidy(dm)
  } else {
    sample_size <- sample_size * get_samples_per_ind(dm)
  }

  sample_size
}




#' @describeIn get_feature_table Returns a vector of groups in the model
#' @export
get_groups <- function(dm) {
  groups <- read_cache(dm, 'groups')
  if (is.null(groups)) {
    groups <- sort(unique(c(sapply(dm$features, function(f) f$get_group()),
                            sapply(dm$sum_stats, function(s) s$get_group()),
                            sapply(dm$loci, function(l) l$get_group()))))

    if (all(groups == 0)) return(0)
    if (!1 %in% groups) groups <- c(1, groups)
    groups <- groups[groups != 0]
    cache(dm, 'groups', groups)
  }
  groups
}


#' @describeIn get_feature_table Returns a matrix with detailed length
#' information about the loci in the model.
#' @export
get_locus_length_matrix <- function(dm, individual=TRUE) {
  llm <- read_cache(dm, "llm")
  if (is.null(llm)) {
    # Select the relevant loci in dm$loci
    loci <- get_loci(dm, idx=TRUE)

    # Repeat the row if number > 1
    if (individual && length(loci) == 1) {
      loci <- rep(loci, dm$loci[[loci]]$get_number())
    }

    assert_that(length(loci) > 0)
    llm <- do.call(rbind, lapply(dm$loci[loci], function(l) l$get_length(TRUE)))
    cache(dm, "llm", llm)
  }
  llm
}


#' @describeIn get_feature_table Returns the population that is marked as outgroup
#' @export
get_outgroup <- function(model) {
  as.integer(search_feature(model, 'outgroup')$parameter)
}


#' @describeIn get_feature_table Returns the number of samples in the outgroup
#' @export
get_outgroup_size <- function(model, for_sim = FALSE) {
  outgroup_size <- get_sample_size(model, for_sim)[get_outgroup(model)]
  if (length(outgroup_size) == 0) outgroup_size <- 0
  outgroup_size
}


#' @describeIn get_feature_table Returns the index of the individuals of one
#'   population
#' @export
get_population_indiviuals <- function(model, pop) {
  if (!pop %in% get_populations(model)) stop('Invalid population')
  sample_size <- get_sample_size(model)
  from <- cumsum(c(0, sample_size)) + 1
  to <- cumsum(sample_size)
  from[pop]:to[pop]
}


get_par_names <- function(model) {
  if (length(get_parameter(model)) == 0) return(character(0))
  sapply(get_parameter(model), function(par) par$get_name())
}


get_loci <- function(model, group=NA, idx=FALSE) {
  if (is.na(group)) loci <- seq(along = model$loci)
  else {
    loci <- which(sapply(model$loci, function(l) l$get_group() == group))
    if (length(loci) == 0) {
      loci <- which(sapply(model$loci, function(l) l$get_group() == 0))
    }
  }
  if (idx) return(loci)
  model$loci[loci]
}
