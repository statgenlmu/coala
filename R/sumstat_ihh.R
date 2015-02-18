#' @importFrom R6 R6Class
#' @importFrom rehh data2haplohh calc_ehh
SumStat_iHH <- R6Class('SumStat_iHH', inherit = SumStat,
  private = list(position = NA),
  public = list(
    calculate = function(seg_sites, files, model) {
      assert_that(is.list(seg_sites))
      assert_that(is.model(model))
      pos <- get_snp_positions(seg_sites, model, relative = FALSE)
      lapply(1:length(seg_sites), function(locus) {
        assert_that(is.matrix(seg_sites[[locus]]))
        ehh <- sapply(1:nrow(seg_sites[[locus]]), function(snp) {
          calc_ehh(self$segsites_to_rehh_data(seg_sites[[locus]], pos[[locus]]),
                   mrk = snp, plotehh = FALSE)$ihh
        })
        colnames(ehh) <- pos[[locus]]
        ehh
      })
    },
    segsites_to_snp_map = function(seg_sites, pos) {
      map <- data.frame(name=seq(along = pos), chr=1, pos = pos, anc=0, der=1)
      file <- tempfile('snp_map')
      write.table(map, file, row.names = FALSE, col.names = FALSE)
      file
    },
    segsites_to_haplo = function(seg_sites) {
      file <- tempfile('haplotypes')
      write.table(cbind(1:ncol(seg_sites), t(seg_sites)), file,
                  row.names = FALSE, col.names = FALSE)
      file
    },
    segsites_to_rehh_data = function(seg_sites, pos) {
      haplo <- self$segsites_to_haplo(seg_sites)
      snp_map <- self$segsites_to_snp_map(seg_sites, pos)
      capture.output(rehh <- data2haplohh(haplo, snp_map, recode.allele = TRUE))
      unlink(c(snp_map, haplo))
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
#' coalsimr relies on the package rehh to calculate this statistic. Please refer
#' to their documentation for detailed information on the concrete
#' implementation.
#'
#' @inheritParams sumstat_file
#' @export
sumstat_iHH <- function(name = 'iHH', group = 0) {
  SumStat_iHH$new(name, group = group)
}
