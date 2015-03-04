sumstat <- R6Class('sumstat',
  private = list(name=NA, group=NA),
  public = list(
    initialize = function(name, group=0) {
      private$name <- name
      private$group <- group
    },
    calculate = function(seg_sites, files, model) {
      stop("Overwrite this function with the calculation of the statistic.")
    },
    get_name = function() private$name,
    get_group = function() private$group
  )
)


is.sum_stat <- function(sum_stat) 'sumstat' %in% class(sum_stat)


# Written to `model$sum_stats` on initialization
create_sumstat_container <- function() {
  list()
}


# Add a summary statistic to a model
add_to_model.sumstat <- function(sum_stat, model, feat_name) {
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
get_summary_statistics <- function(dm, group = 1) {
  ss <- sapply(get_group_statistics(dm, group), function(stat) stat$get_name())
  if (length(ss) == 0) return(character())
  names(ss) <- NULL
  ss
}


get_group_statistics <- function(model, group) {
  if (group == 'all') return(model$sum_stats)
  stats <- lapply(model$sum_stats, function(sum_stat) {
    if (sum_stat$get_group() %in% c(0, group)) return(sum_stat)
    else return(NA)
  })
  stats[!is.na(stats)]
}


get_sumstat_groups <- function(model) {
  unique(sapply(model$sum_stats, function(ss) ss$get_group()))
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
