sumstat_sgtrees_class <- R6Class('sumstat_Trees', inherit = sumstat,
  private = list(loci_length = NULL),
  public = list(
    initialize = function(loci_length, group) {
      private$loci_length <- loci_length
      super$initialize('sg_trees', group)
    },
    calculate = function(seg_sites, files, model) {
      parse_trees(files[[1]], private$loci_length, tempfile)
    }
  )
)

#' Returns ancestral tress in NEWICK fromat from simulations
#'
#' @inheritParams sumstat_file
sumstat_sg_trees <- function(locus_length, group = 0) {
  Feature$new('trees', par_const(NA), group = group) +
    sumstat_sgtrees_class$new(locus_length, group)
}
