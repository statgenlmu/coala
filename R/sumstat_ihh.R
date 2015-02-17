#' @importFrom R6 R6Class
#' @importFrom rehh data2haplohh calc_ehh
SumStat_iHH <- R6Class('SumStat_iHH', inherit = SumStat,
  public = list(
    calculate = function(seg_sites, files, model) {
      calc_ehh(self$segsites_to_rehh_data(seg_sites, model),
               mrk = 2, plotehh = FALSE)
    },
    segsites_to_snp_map = function(seg_sites, model) {
      assert_that(is.matrix(seg_sites))
      if (has_trios(model, private$group))
        stop('iHH does not support trios yet.')
      pos <- attr(seg_sites, 'positions') * get_locus_length(model,
                                                             private$group)
      map <- data.frame(name=seq(along = pos), chr=1, pos = pos, anc=0, der=1)

      file <- tempfile('snp_map')
      write.table(map, file, row.names = FALSE, col.names = FALSE)
      file
    },
    segsites_to_haplo = function(seg_sites) {
      assert_that(is.matrix(seg_sites))
      file <- tempfile('haplotypes')
      write.table(cbind(1:ncol(seg_sites), t(seg_sites)), file,
                  row.names = FALSE, col.names = FALSE)
      file
    },
    segsites_to_rehh_data = function(seg_sites, model) {
      output <- tempfile('rehh_output')
      sink(output)
      haplo <- self$segsites_to_haplo(seg_sites)
      snp_map <- self$segsites_to_snp_map(seg_sites, model)
      rehh <- data2haplohh(haplo, snp_map, recode.allele = TRUE)
      sink(NULL)
      unlink(c(haplo, snp_map, output))
      rehh
    }
  )
)

#' Integrated Extended Haplotype Homozygosity
#'
#' This summary statistic calculates a modified version of the iHH statistic
#' introduced by
#'
#'  Voight et al., A map of recent positive selection in the human genome.
#'  PLoS Biol, 4(3):e72, Mar 2006
#'
#' their documentation for detailed information.
#'
#' @inheritParams sumstat_file
#' @export
sumstat_iHH <- function(name = 'iHH', group = 0) {
  SumStat_iHH$new(name, group = group)
}
