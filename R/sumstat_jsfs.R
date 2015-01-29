#' Returns the Joint Site Frequency Spectrum from simulations
#'
#' @inheritParams sumstat_file
#' @export
sumstat_jsfs <- function(group = 0) {
  SumStat$new('jsfs', group = group)
}
