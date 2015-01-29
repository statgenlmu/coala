#' Returns ancestral tress in NEWICK fromat from simulations
#'
#' @inheritParams sumstat_file
#' @export
sumstat_trees <- function(group = 0) {
  SumStat$new('trees', group = group)
}
