#' @importFrom R6 R6Class
sumstat_dna_class <- R6Class('sumstat_DNA', inherit = sumstat,
  public = list(
    calculate = function(seg_sites, files, model) seg_sites,
    parse_files = function(files) {
      lapply(files, function(file) {
        scan(file, what="character")
      })
    }
  )
)

#' Returns the Segregation Sites Statistics from simulations
#'
#' @inheritParams sumstat_file
#' @export
sumstat_seg_sites <- function(name = 'dna', group = 0) {
  sumstat_dna_class$new(name, group = group)
}
