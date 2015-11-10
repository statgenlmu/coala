#' Base Class for Summary Statistics
#'
#' If you want to create additional summary statistics for coala, create
#' \pkg{R6} classes that inherit from this object.
#'
#' @export
sumstat_class <- R6Class("sumstat",
  private = list(
    name = NA,
    req_files = FALSE,
    req_trees = FALSE,
    req_segsites = FALSE,
    transformation = NULL
  ),
  public = list(
    initialize = function(name, transformation) {
      assert_that(is.character(name))
      assert_that(length(name) == 1)
      private$name <- name

      assert_that(is.function(transformation))
      private$transformation <- transformation
    },
    calculate = function(seg_sites, trees, files, model) {
      stop("Overwrite this function with the calculation of the statistic.")
    },
    check = function(model) {
      # Optional functions that checks if a model is compatible with the stat.
      # Should throw an informative error if not.
      # This gets executed before the stat is added to the model.
      invisible(TRUE)
    },
    get_name = function() private$name,
    requires_files = function() private$req_files,
    requires_segsites = function() private$req_segsites,
    requires_trees = function() private$req_trees,
    print = function() cat(class(self)[1], "\n"),
    transform = function(x) private$transformation(x)
  )
)


is.sum_stat <- function(sumstat) inherits(sumstat, "sumstat")


# Written to `model$sum_stats` on initialization
create_sumstat_container <- function() list()


# Add a summary statistic to a model
add_to_model.sumstat <- function(sum_stat, model, feat_name) {
  sum_stat$check(model)

  if (sum_stat$get_name() %in% names(model$sum_stats))
    stop("Can not add ", feat_name, " to model: ",
         "There is already a statistic with name ", sum_stat$get_name())

  if (sum_stat$requires_files() && !requires_files(model))
    model <- model + files_feat_class$new()
  if (sum_stat$requires_segsites() && !requires_segsites(model))
    model <- model + segsites_feat_class$new()
  if (sum_stat$requires_trees() && !requires_trees(model))
    model <- model + trees_feat_class$new()

  # Save the statistic
  model$sum_stats[[sum_stat$get_name()]] <- sum_stat

  # Update cache
  model$id <- get_id()
  model
}


#' @param pop The population for which aspects are returned
#' @describeIn get_features Returns the summary statistics in the model
#' @export
get_summary_statistics <- function(model) {
  model$sum_stats
}


calc_sumstats <- function(seg_sites, trees, files, model,
                          pars, cmds, simulator) {
  if (missing(pars)) pars <- numeric(0)
  stopifnot(is.model(model))

  if (is.list(cmds) && simulator$get_name() != "seqgen") {
    cmds <- do.call(c, cmds)
  }

  sum_stats <- list(pars = pars,
                    cmds = cmds,
                    simulator = simulator$get_info())

  # Process seg_sites for trios and unphase if neccessary
  if (requires_segsites(model)) {
    if (has_trios(model) && simulator$get_name() != "seqgen") {
      seg_sites <- conv_for_trios(seg_sites, model)
    }

    if (is_unphased(model)) {
      seg_sites <- unphase_segsites(seg_sites,
                                    get_ploidy(model),
                                    get_samples_per_ind(model))
    }

    assert_that(all(vapply(seg_sites, is_segsites, logical(1))))
  }

  for (stat in model$sum_stats) {
    sum_stats[[stat$get_name()]] <-
      stat$transform(stat$calculate(seg_sites, trees, files, model))
  }

  sum_stats
}

requires_segsites <- function(model) {
  any(sapply(get_summary_statistics(model), function(x) x$requires_segsites()))
}

requires_trees <- function(model) {
  any(sapply(get_summary_statistics(model), function(x) x$requires_trees()))
}

requires_files <- function(model) {
  any(sapply(get_summary_statistics(model), function(x) x$requires_files()))
}
