#' @importFrom R6 R6Class
SumstatIhh <- R6Class('sumstat_ihh', inherit = Sumstat, #nolint
  private = list(
    req_segsites = TRUE,
    position = NA,
    population = NULL,
    get_snp = function(positions, locus, model) {
      if (is.na(private$position)) return(seq(along = positions))
      pos <- conv_middle_to_trio_pos(private$position, model,
                                     relative_out = FALSE)[locus]
      which.min(abs(pos - positions))
    }),
  public = list(
    initialize = function(name, position, population) {
      if (!requireNamespace("rehh", quietly = TRUE)) {
        stop("Package rehh is required to calculate the iHH summary statistic.",
             " Please install it.", call. = FALSE)
      }
      assert_that(is.numeric(population))
      assert_that(length(population) == 1)
      private$population <- population
      private$position <- position
      super$initialize(name)
    },
    calculate = function(seg_sites, files, model) {
      assert_that(is.list(seg_sites))
      assert_that(is.model(model))
      pos <- get_snp_positions(seg_sites, model, relative = FALSE)
      ind <- get_population_indiviuals(model, private$population)
      lapply(1:length(seg_sites), function(locus) {
        assert_that(is.matrix(seg_sites[[locus]]))
        if (ncol(seg_sites[[locus]]) == 0) {
          return(matrix(0, 2, 0, dimnames = list(c("Anc. Allele",
                                                   "Der. Allele"), NULL)))
        }
        snps <- private$get_snp(pos[[locus]], locus, model)
        ehh <- sapply(snps, function(snp) {
          rehh::calc_ehh(self$segsites_to_rehh_data(seg_sites[[locus]],
                                                    pos[[locus]],
                                                    ind),
                         mrk = snp, plotehh = FALSE)$ihh
        })
        colnames(ehh) <- pos[[locus]][snps]
        ehh
      })
    },
    segsites_to_snp_map = function(seg_sites, pos) {
      map <- data.frame(name=seq(along = pos), chr=1, pos = pos, anc=0, der=1)
      file <- tempfile('snp_map')
      write.table(map, file, row.names = FALSE, col.names = FALSE)
      file
    },
    segsites_to_haplo = function(seg_sites, ind) {
      file <- tempfile('haplotypes')
      write.table(cbind(ind, seg_sites[ind, ]), file,
                  row.names = FALSE, col.names = FALSE)
      file
    },
    segsites_to_rehh_data = function(seg_sites, pos, ind) {
      haplo <- self$segsites_to_haplo(seg_sites, ind)
      snp_map <- self$segsites_to_snp_map(seg_sites, pos)
      capture.output(rehh <- rehh::data2haplohh(haplo,
                                                snp_map,
                                                recode.allele = TRUE))
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
#' implementation. It is required to install the package \code{rehh} to use this
#' function.
#'
#' @inheritParams sumstat_file
#' @param position A position relative to the locus extent, e.g. 0.5 for the
#'   middle of the locus. If provided, the iHH will be calculate using the SNP
#'   closest to the given position as focal SNP. Otherwise, all SNPs will be
#'   used as focal SNPs in turn, and all the values reportet. If trios are used,
#'   than the position is relative to the middle locus' extend.
#' @export
sumstat_ihh <- function(name = 'ihh', position=NA, population=1) {
  SumstatIhh$new(name, position, population) #nolint
}
