#' @importFrom R6 R6Class
SumstatSegSites <- R6Class('SumstatSegSites', inherit = Sumstat, #nolint
  private = list(req_segsites = TRUE),
  public = list(
    calculate = function(seg_sites, files, model) seg_sites
  )
)

#' Returns the Segregation Sites Statistics from simulations
#'
#' @inheritParams sumstat_file
#' @export
sumstat_seg_sites <- function(name = 'seg_sites') {
  SumstatSegSites$new(name) #nolint
}
