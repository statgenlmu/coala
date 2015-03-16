#' @importFrom R6 R6Class
SumstatDna <- R6Class('SumstatDna', inherit = Sumstat, #nolint
  private = list(req_files = TRUE),
  public = list(
    calculate = function(seg_sites, files, model) {
      dna <- parse_sg_output(files,
                            sum(get_sample_size(model, for_sim = TRUE)),
                            get_locus_length_matrix(model),
                            get_locus_number(model),
                            calc_seg_sites = FALSE)
      lapply(dna, function(locus) {
        apply(locus, 2, function(x) attr(locus, 'levels')[x])
      })
    }
  )
)

#' Returns the Segregation Sites Statistics from simulations
#'
#' @inheritParams sumstat_file
#' @export
sumstat_dna <- function(name = 'dna') {
  coal_model() +
    Feature$new('sumstat_dna', par_const(NA)) +
    SumstatDna$new(name) #nolint
}
