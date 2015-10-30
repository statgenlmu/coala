#' @importFrom R6 R6Class
stat_dna_class <- R6Class("dna_stat", inherit = sumstat_class,
  private = list(),
  public = list(
    calculate = function(seg_sites, trees, dna, model) {
      if (requires_segsites(model)) {
        stop("Can not generate both seg. sites and DNA")
      }
      if (has_trios(model)) {
        stop("Cannot generate DNA for trio models")
      }
      lapply(dna, function(locus) {
        apply(locus, 2, function(x) attr(locus, "levels")[x])
      })
    }
  )
)

#' Returns the Segregation Sites Statistics from simulations
#'
#' @inheritParams sumstat_four_gamete
#' @export
sumstat_dna <- function(name = "dna", transformation = identity) {
  stat_dna_class$new(name, transformation)
}
