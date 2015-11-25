ign_singletons_class <- R6Class("ign_singletons", inherit = feature_class)


#' Feature: Ignore Singletons
#'
#' Mutations that are observed in just one haplotype ('singletons') are often
#' regarded as likely candidates for sequencing errors. Sometimes, it can be
#' advantageous to exclude them from an analysis. This feature removes all
#' singletons from the simulated data before the summary statistics are
#' calculated.
#'
#' For this function, a singleton a mutation for which the derived allele is
#' observed exactly once in all populations.
#'
#' @return The feature, which can be added to a model using `+`.
#' @export
#' @examples
#' model <- coal_model(2, 1) +
#'   feat_mutation(10) +
#'   feat_ignore_singletons() +
#'   sumstat_sfs("n_mut", transformation = sum)
#' simulate(model)$n_mut # Is 0, because all SNPs are singletons in this model
feat_ignore_singletons <- function() {
  ign_singletons_class$new()
}


has_ign_singletons <- function(model) {
  any(vapply(get_features(model), inherits,
             what = "ign_singletons", logical(1)))
}


remove_singletons <- function(segsites_list) {
  lapply(segsites_list, function(segsites) {
    is_singleton <- colSums(get_snps(segsites)) == 1
    segsites[ , !is_singleton]
  })
}


#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_ms_arg.ign_singletons <- ignore_par

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_msms_arg.ign_singletons <- ignore_par

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_scrm_arg.ign_singletons <- ignore_par

#' @describeIn conv_to_ms_arg Feature conversion
#' @export
conv_to_seqgen_arg.ign_singletons <- ignore_par
