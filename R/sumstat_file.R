#' Returns files with the raw results of simulations
#'
#' @param group The locus group for which this summary statistic is reported.
#'   The default of `0` corresponds to all groups.
#' @export
sumstat_file <- function(group = 0) {
  SumStat$new('file', group = group)
}
