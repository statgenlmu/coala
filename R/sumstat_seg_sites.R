#' @importFrom R6 R6Class
sumstat_segsites_class <- R6Class('sumstat_SegSites', inherit = sumstat,
  public = list(
    calculate = function(seg_sites, files, model) seg_sites
  )
)

#' Returns the Segregation Sites Statistics from simulations
#'
#' @inheritParams sumstat_file
#' @export
sumstat_seg_sites <- function(name = 'seg_sites', group = 0) {
  sumstat_segsites_class$new(name, group = group)
}
