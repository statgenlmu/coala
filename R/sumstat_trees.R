sumstat_sgtrees_class <- R6Class('sumstat_Trees', inherit = sumstat,
  public = list(
    calculate = function(seg_sites, files, model) {
      llm <- get_locus_length_matrix(model, has_inter_locus_var(model))
      assert_that(length(files) == nrow(llm))
      lapply(1:nrow(llm), function(i) {
        parse_trees(files[[i]], llm[i, 1:5], tempfile)
      })
    }
  )
)

#' Returns ancestral tress in NEWICK fromat from simulations
#'
#' @inheritParams sumstat_file
sumstat_sg_trees <- function() {
  coal_model() +
    Feature$new('trees', par_const(NA)) +
    sumstat_sgtrees_class$new('trees')
}
