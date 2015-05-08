Feature_segsites <- R6Class("Feature_seg_sites", inherit = Feature,
  public = list(
    print = function() cat("Generating Seg. Sites\n")
  )
)
conv_to_ms_arg.Feature_seg_sites <- function(feature, model) {
  if (!any(vapply(get_features(model), is_feat_mutation, logical(1)))) {
    stop("model requires mutation to calculate summary statistics")
  }
  ""
}
conv_to_msms_arg.Feature_seg_sites <- conv_to_ms_arg.Feature_seg_sites
conv_to_scrm_arg.Feature_seg_sites <- conv_to_ms_arg.Feature_seg_sites
conv_to_seqgen_arg.Feature_seg_sites <- conv_to_ms_arg.Feature_seg_sites

Feature_trees <- R6Class("Feature_trees", inherit = Feature,
  public = list(
    print = function() cat("Generating Trees\n")
  )
)
conv_to_ms_arg.Feature_trees <- function(feature, model) "-T"
conv_to_msms_arg.Feature_trees <- conv_to_ms_arg.Feature_trees
conv_to_scrm_arg.Feature_trees <- conv_to_ms_arg.Feature_trees
conv_to_seqgen_arg.Feature_trees <- function(feature, model) {
  stop("generation of trees is not supported.")
}

Feature_files <- R6Class("Feature_files", inherit = Feature,
  public = list(
    print = function() cat("Generating Files\n")
  )
)
conv_to_ms_arg.Feature_files <- ignore_par
conv_to_msms_arg.Feature_files <- ignore_par
conv_to_scrm_arg.Feature_files <- ignore_par
conv_to_seqgen_arg.Feature_files <- ignore_par
