stat_trees_class <- R6Class("stat_trees", inherit = sumstat_class,
  private = list(req_files = FALSE, req_trees = TRUE),
  public = list(
    calculate = function(seg_sites, trees, files, model, sim_tasks = NULL) {
      trees
    }
  )
)

#' Summary Statistic: Ancestral Trees
#'
#' This statistic returns ancestral tress in NEWICK format.
#'
#' @export
#' @inheritParams sumstat_four_gamete
#' @template summary_statistics
#' @examples
#' # Without recombination:
#' model <- coal_model(4, 2) + sumstat_trees()
#' stats <- simulate(model)
#' print(stats$trees)
#'
#' # With recombination:
#' model <- model + feat_recombination(5)
#' stats <- simulate(model)
#' print(stats$trees)
sumstat_trees <- function(name = "trees") {
  stat_trees_class$new(name, identity)
}
