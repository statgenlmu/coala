#' @importFrom R6 R6Class
stat_jsfs_class <- R6Class("stat_jsfs", inherit = sumstat_class,
  private = list(
    populations = NULL,
    req_segsites = TRUE,
    per_locus = FALSE
  ),
  public = list(
    initialize = function(name, populations, per_locus, transformation) {
      assert_that(is.numeric(populations))
      private$populations <- populations

      assert_that(is.logical(per_locus))
      assert_that(length(per_locus) == 1)
      private$per_locus <- per_locus
      super$initialize(name, transformation)
    },
    calculate = function(seg_sites, trees, files, model) {
      ind_per_pop <- lapply(private$populations, get_population_indiviuals,
                            model = model)

      if (private$per_locus) {
        jsfs <- lapply(seg_sites, function(x) {
          calc_jsfs(list(x), ind_per_pop)
        })
      } else {
        jsfs <- calc_jsfs(seg_sites, ind_per_pop)
      }
      jsfs
    }
  )
)

#' Calculates the Joint Site Frequency Spectrum from simulations
#'
#' @inheritParams sumstat_four_gamete
#' @param populations The populations for which the statistic is calculated.
#' @param per_locus If \code{TRUE}, the JSFS is return for each locus instead
#'   of globally. In this case, the result is a list, where each entry is the
#'   JSFS for the corresponding locus.
#' @export
#' @examples
#' set.seed(75)
#' model <- coal_model(2:4, 2) +
#'   feat_mutation(5) +
#'   feat_migration(1, symmetric = TRUE) +
#'   sumstat_jsfs("jsfs_12", populations = c(1, 2)) +
#'   sumstat_jsfs("jsfs_123", populations = 1:3)
#' stats <- simulate(model)
#' print(stats$jsfs_12)
#' print(stats$jsfs_123)
sumstat_jsfs <- function(name = "jsfs", populations = c(1, 2),
                         per_locus = FALSE, transformation = identity) {
  stat_jsfs_class$new(name, populations, per_locus, transformation)
}
