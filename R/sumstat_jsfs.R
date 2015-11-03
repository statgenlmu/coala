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
      assert_that(length(populations) == 2)
      private$populations <- populations

      assert_that(is.logical(per_locus))
      assert_that(length(per_locus) == 1)
      private$per_locus <- per_locus
      super$initialize(name, transformation)
    },
    calculate = function(seg_sites, trees, files, model) {
      ind_pop1 <- get_population_indiviuals(model, private$populations[1])
      ind_pop2 <- get_population_indiviuals(model, private$populations[2])

      if (private$per_locus) {
        jsfs <- lapply(seg_sites, function(x) {
          calc_jsfs(list(x), ind_pop1, ind_pop2)
        })
      } else {
        jsfs <- calc_jsfs(seg_sites, ind_pop1, ind_pop2)
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
sumstat_jsfs <- function(name = "jsfs", populations = c(1, 2),
                         per_locus = FALSE, transformation = identity) {
  stat_jsfs_class$new(name, populations, per_locus, transformation)
}
