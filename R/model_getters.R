#' Getters for coalescent models
#'
#' @param model The coalescent model from which aspects are returned
#'
#' @export
#' @author Paul Staab
#' @describeIn get_feature_table Returns the features of a model as a data.frame
get_feature_table <- function(model) {
  stopifnot(is.model(model))
  stopifnot(!is.null(model$features))
  model$features
}


#' @export
#' @describeIn get_feature_table Returns the ranged parameters of a model as a
#'   data.frame
get_parameter_table <- function(model) {
  stopifnot(is.model(model))
  stopifnot(!is.null(model$parameters))
  model$parameters
}


#' @param group The locus group for which aspects are returned
#' @export
#' @describeIn get_feature_table Returns the length of the loci in a locus group
get_locus_length <- function(dm, group=1) {
  length <- dm$loci[dm$loci$group == group, c(6,8,10), drop = FALSE]
  if (nrow(length) == 0) {
    length <- dm$loci[dm$loci$group == 0, c(6,8,10), drop = FALSE]
  }
  if (nrow(length) == 0) stop("Failed to determine loci length")
  as.integer(rowSums(length))
}


#' @describeIn get_feature_table Returns the number of loci in a locus group
#' @export
get_locus_number <- function(dm, group=1) {
  number <- dm$loci[dm$loci$group == group, 'number']
  if (length(number) == 0) {
    if (group > 0) number <- dm$loci[dm$loci$group == 0, 'number']
    else stop("Failed to determine loci number")
  }
  as.integer(sum(number))
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
  if (all(c(dm$features$group == 0,
            dm$sum_stats$group == 0,
            dm$loci$group == 0))) return(1)

  groups <- sort(unique(c(1,
                          dm$features$group,
                          dm$sum_stats$group,
                          dm$loci$group)))

  groups[groups != 0]
}


#' @describeIn get_feature_table Returns a matrix with detailed length
#' information about the loci in the model.
#' @export
get_locus_length_matrix <- function(dm, group=1) {
  # Select the rows of the group
  rows <- which(dm$loci$group == group)
  if (sum(rows) == 0) rows <- which(dm$loci$group == 0)

  # Repeat the row if number > 1
  if (length(rows) == 1) rows <- rep(rows, dm$loci$number[rows])

  # Return the matrix
  llm <- dm$loci[rows, 6:10, drop = FALSE]
  row.names(llm) <- NULL
  as.matrix(llm)
}


#' @describeIn get_feature_table Returns the population that is marked as outgroup
#' @export
get_outgroup <- function(model) {
  as.integer(search_feature(model, 'outgroup')$parameter)
}


#' @describeIn get_feature_table Returns the number of samples in the outgroup
#' @export
get_outgroup_size <- function(model, for_sim = FALSE) {
  get_sample_size(model, for_sim)[get_outgroup(model)]
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
