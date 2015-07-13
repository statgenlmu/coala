#' @importFrom R6 R6Class
stat_dna_class <- R6Class("dna_stat", inherit = sumstat_class,
  private = list(req_files = TRUE),
  public = list(
    calculate = function(seg_sites, trees, files, model) {
      dna <- parse_sg_output(files,
                            sum(get_sample_size(model, for_sim = TRUE)),
                            get_locus_length_matrix(model),
                            get_locus_number(model),
                            calc_seg_sites = FALSE)
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
sumstat_dna <- function(name = "dna") {
  stat_dna_class$new(name) #nolint
}
