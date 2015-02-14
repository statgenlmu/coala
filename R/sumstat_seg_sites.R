#' @importFrom R6 R6Class
SumStat_SegSites <- R6Class('SumStat_SegSites', inherit = SumStat,
  public = list(
    calculate = function(seg_sites, files, model) seg_sites
  )
)

#' Returns the Segregation Sites Statistics from simulations
#'
#' @inheritParams sumstat_file
#' @export
sumstat_seg_sites <- function(name = 'seg_sites', group = 0) {
  SumStat_SegSites$new(name, group = group)
}
