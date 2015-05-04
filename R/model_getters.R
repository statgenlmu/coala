#' Getters for coalescent models
#'
#' @param model The coalescent model from which aspects are returned
#'
#' @export
#' @author Paul Staab
get_features <- function(model) model$features


#' @export
#' @describeIn get_features Returns the ranged parameters of a model as a
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
#' @describeIn get_features Returns the ranged parameters of a model
get_parameter <- function(model) {
  stopifnot(is.model(model))
  model$parameter
}


#' @param locus The number of the locus.
#' @param total If \code{FALSE}, the length of loci in a trio will be reported
#'   individually. If \code{TRUE} the sum of the loci's length will be reported.
#'   This does not affect non-trio loci.
#'
#' @describeIn get_features Returns the length of the loci in a locus group
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


#' @describeIn get_features Returns a vector of populations in the model
#' @export
get_populations <- function(model) {
  unique(search_feature(model, 'sample')$pop.source)
}


#' @describeIn get_features Returns a matrix with detailed length
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


#' @describeIn get_features Returns the number of loci in a locus group
#' @export
get_locus_number <- function(model, group=NA) {
  numbers <- get_locus_length_matrix(model)[ , "number"]
  if (is.na(group)) return(sum(numbers))
  if (has_inter_locus_var(model)) return(1)
  numbers[group]
}


#' @describeIn get_features Returns the population that is marked as outgroup
#' @export
get_outgroup <- function(model) {
  as.integer(search_feature(model, 'outgroup')$parameter)
}


#' @describeIn get_features Returns the number of samples in the outgroup
#' @export
get_outgroup_size <- function(model, for_sim = FALSE) {
  outgroup_size <- get_sample_size(model, for_sim)[get_outgroup(model)]
  if (length(outgroup_size) == 0) outgroup_size <- 0
  outgroup_size
}


#' @describeIn get_features Returns the index of the individuals of one
#'   population
#' @export
get_population_indiviuals <- function(model, pop, zero_indexed = FALSE) {
  if (pop == "all") return(1:sum(get_sample_size(model)))

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
