Feature_unphased <- R6Class("Feature_unphased", inherit = Feature,
  private = list(ploidy = NA, samples_per_ind = NA),
  public = list(
    initialize = function(ploidy, samples_per_ind) {
      assert_that(length(ploidy) == 1)
      assert_that(is.numeric(ploidy))
      assert_that(length(samples_per_ind) == 1)
      assert_that(is.numeric(samples_per_ind))
      assert_that(samples_per_ind <= ploidy)
      private$ploidy = ploidy
      private$samples_per_ind = samples_per_ind
    },
    get_ploidy = function() private$ploidy,
    get_samples_per_ind = function() private$samples_per_ind,
    print = function() {
      cat("Unphasing of", private$ploidy, "chromosomes into",
          private$samples_per_ind, "pseudo-chromosomes\n")
    }
  )
)

#' Feature: Unphased
#'
#' This simulated unphased data by creating "pseudo-chromosomes". For these,
#' each position is randomly taken from a phased chromosome obtained by
#' simulation.
#'
#' If this is used, the sample size is understood as the number of individuals.
#' For each individual, \code{poidy} chromosomes are simulated, and
#' \code{samples_per_ind} pseudo-chromosomes are created of these.
#'
#' @param ploidy The number of phased chromosomes that are simulated per
#'   individual.
#' @param samples_per_ind The number of pseudo-chromosomes that are created
#'   from the phased chromosomes for each indidivual.
#' @return The feature, which can be added to a model using `+`.
#' @export
feat_unphased <- function(ploidy, samples_per_ind=ploidy) {
  Feature_unphased$new(ploidy, samples_per_ind)
}


is_feat_unphased <- function(feat) any("Feature_unphased" == class(feat))

get_feature_unphased <- function(model) {
  mask <- vapply(model$features, is_feat_unphased, logical(1))
  if (!any(mask)) return(NULL)
  if (sum(mask) > 1) stop("multiple unphasings are not supported")
  model$features[mask][[1]]
}


get_ploidy <- function(model) {
  feat <- get_feature_unphased(model)
  if (is.null(feat)) return(1L)
  as.integer(feat$get_ploidy())
}


get_samples_per_ind <- function(model) {
  feat <- get_feature_unphased(model)
  if (is.null(feat)) return(1L)
  as.integer(feat$get_samples_per_ind())
}


is_unphased <- function(model) !is.null(get_feature_unphased(model))

conv_to_ms_arg.Feature_unphased <- ignore_par
conv_to_msms_arg.Feature_unphased <- ignore_par
conv_to_seqgen_arg.Feature_unphased <- ignore_par
