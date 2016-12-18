#' @importFrom R6 R6Class
stat_dna_class <- R6Class("dna_stat", inherit = sumstat_class,
  private = list(req_files = TRUE),
  public = list(
    calculate = function(seg_sites, trees, files, model, sim_tasks) {
      if (has_trios(model)) {
        stop("Cannot generate DNA for trio models")
      }
      if (is_unphased(model)) {
        stop("Cannot generate unphased DNA")
      }

      assert_that(length(files) == length(sim_tasks))

      lapply(seq_along(files), function(i) {
        assert_that(file.exists(files[i]))
        locus_length <- sum(sim_tasks[[i]]$get_arg("trio_dists"))
        dna <- parse_seqgen_output(readLines(files[i]),
                                   individuals = sum(get_sample_size(model,
                                                                     TRUE)),
                                   locus_length = locus_length,
                                   locus_number = sim_tasks[[i]]$locus_number,
                                   outgroup_size = get_outgroup_size(model,
                                                                     TRUE),
                                   calc_segsites = FALSE)
        dna[[1]]
      })
    }
  )
)

#' Summary Statistic: DNA
#'
#' This summary statistic returns the actual DNA sequences from
#' finite sites simulations. It can not be
#' calculated together with other summary statistics or when assuming
#' an infinite sites mutation model. No outgroup
#' is needed for it, and the outgroup sequences will also be
#' returned if present.
#'
#' @inheritParams sumstat_four_gamete
#' @return A list of sequences for each locus. Each entries is a
#'         character matrix decoding the sequences. Each row
#'         is an individual, and each column is a genetic position.
#' @template summary_statistics
#' @export
#' @examples
#' model <- coal_model(5, 1, 10) +
#'  feat_mutation(5, model = "GTR", gtr_rates = rep(1, 6)) +
#'  sumstat_dna()
#' \dontrun{simulate(model)$dna}
sumstat_dna <- function(name = "dna", transformation = identity) {
  stat_dna_class$new(name, transformation)
}
