Sumstat <- R6Class('Sumstat',
  private = list(
    name=NA,
    req_files=FALSE,
    req_segsites=FALSE
  ),
  public = list(
    initialize = function(name) {
      private$name <- name
    },
    calculate = function(seg_sites, files, model) {
      stop("Overwrite this function with the calculation of the statistic.")
    },
    get_name = function() private$name,
    get_group = function() private$group,
    requires_files = function() private$req_files,
    requires_segsites = function() private$req_segsites
  )
)


is.sum_stat <- function(sum_stat) 'Sumstat' %in% class(sum_stat)


# Written to `model$sum_stats` on initialization
create_sumstat_container <- function() {
  list()
}


# Add a summary statistic to a model
add_to_model.Sumstat <- function(sum_stat, model, feat_name) {
  if (sum_stat$get_name() %in% names(model$sum_stats))
    stop("Can't add ", feat_name, " to model: ",
         "There is already a statistic with name ", sum_stat$get_name())

  # Save the statistic
  model$sum_stats[[sum_stat$get_name()]] <- sum_stat
  model$id <- get_id()
  model
}


#' @param pop The population for which aspects are returned
#' @describeIn get_feature_table Returns the summary statistics in the model
#' @export
get_summary_statistics <- function(model) {
  model$sum_stats
}


calc_sumstats <- function(seg_sites, files, model, pars) {
  stopifnot(is.model(model))
  sum_stats <- list()
  if (!missing(pars)) sum_stats[['pars']] <- pars

  if (is_unphased(model)) {
    seg_sites <- unphase_segsites(seg_sites,
                                  get_ploidy(model),
                                  get_samples_per_ind(model))
  }

  for (stat in model$sum_stats) {
    sum_stats[[stat$get_name()]] <- stat$calculate(seg_sites, files, model)
  }

  sum_stats
}

requires_segsites <- function(model) {
  any(sapply(get_summary_statistics(model), function(x) x$requires_segsites()))
}

requires_files <- function(model) {
  any(sapply(get_summary_statistics(model), function(x) x$requires_files()))
}
