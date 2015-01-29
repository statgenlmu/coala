#' Returns the Segregation Sites Statistics from simulations
#'
#' @inheritParams sumstat_file
#' @export
sumstat_seg_sites <- function(group = 0) {
  SumStat$new('seg.sites', group = group)
}
