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
    if (!all(sapply(get_parameter(model), is.ranged_par))) {
      stop("Can not create a parameter table with non-ranged pars in model")
    }
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


#' @param locus The number of the locus.
#' @param total If \code{FALSE}, the length of loci in a trio will be reported
#'   individually. If \code{TRUE} the sum of the loci's length will be reported.
#'   This does not affect non-trio loci.
#'
#' @describeIn get_feature_table Returns the length of the loci in a locus group
#' @export
get_locus_length <- function(model, locus=NULL, group=NULL, total=TRUE) {
  assert_that(!(is.null(locus) & is.null(group)))
  llm <- get_locus_length_matrix(model)

  # Group and locus are identical for ilv models
  if (!is.null(group) && has_inter_locus_var(model)) {
    locus <- group
  }

  # Determine to which group the locus belongs
  if (!is.null(locus)) {
    group <- get_locus_group(model, locus)
  }

  ll <- llm[group, 1:5]
  if (total) return(sum(ll))
  if (sum(ll[-3]) == 0) return(ll[3])
  ll
}


get_locus_group <- function(model, locus) {
  llm <- get_locus_length_matrix(model)
  min(which(cumsum(llm[ , "number"]) >= locus))
}


get_locus_group_number <- function(model) {
  if (has_inter_locus_var(model)) return(get_locus_number(model))
  nrow(get_locus_length_matrix(model))
}


#' @describeIn get_feature_table Returns a vector of populations in the model
#' @export
get_populations <- function(model) {
  unique(search_feature(model, 'sample')$pop.source)
}


#' @describeIn get_feature_table Returns a vector of samples sizes per
#'   population.
#' @param for_sim If true, the sample size used internally for the simulation
#'   will be reported rather than the number of actuall samples. The numbers
#'   can be unequal for the simulation of unphased data.
#' @export
get_sample_size <- function(model, for_sim=FALSE) {
  sample_size <- read_cache(model, paste0("sample_size_", for_sim))

  if (is.null(sample_size)) {
    feat.samples <- search_feature(model, type="sample")
    stopifnot(nrow(feat.samples) > 0)

    sample_size <- rep(0, length(get_populations(model)))
    for (row.nr in 1:nrow(feat.samples)) {
      stopifnot(sample_size[feat.samples$pop.source[row.nr]] == 0)
      sample_size[feat.samples$pop.source[row.nr]] <-
        as.integer(feat.samples$parameter[row.nr])
    }

    if (for_sim) {
      sample_size <- sample_size * get_ploidy(model)
    } else {
      sample_size <- sample_size * get_samples_per_ind(model)
    }
    cache(model, paste0("sample_size_", for_sim), sample_size)
  }

  sample_size
}



#' @describeIn get_feature_table Returns a matrix with detailed length
#' information about the loci in the model.
#' @export
get_locus_length_matrix <- function(model) {
  llm <- read_cache(model, "llm")
  if (is.null(llm)) {
    assert_that(length(model$loci) >= 0)
    llm <- do.call(rbind, lapply(model$loci, function(l) {
             number <- ifelse(l$get_number() > 1,
                              round(l$get_number() / model$scaling_factor),
                              l$get_number())
             c(l$get_length(TRUE), number = number)
           }))

    cache(model, "llm", llm)
  }
  llm
}


#' @describeIn get_feature_table Returns the number of loci in a locus group
#' @export
get_locus_number <- function(model, group=NA) {
  numbers <- get_locus_length_matrix(model)[ , "number"]
  if (is.na(group)) return(sum(numbers))
  if (has_inter_locus_var(model)) return(1)
  numbers[group]
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


get_par_names <- function(model, without_priors=FALSE) {
  param <- get_parameter(model)
  if (length(param) == 0) return(character(0))
  if (without_priors) {
    param <- param[!vapply(param, is.prior_par, numeric(1))]
  }
  sapply(param, function(par) par$get_name())
}

get_loci <- function(model) model$loci

get_cmd <- function(model) select_simprog(model)$get_cmd(model)
