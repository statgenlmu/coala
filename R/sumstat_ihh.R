#' @importFrom R6 R6Class
#' @importClassesFrom rehh haplohh
#' @importFrom rehh calc_ehh calc_ehhs
#' @importFrom methods new checkAtAssignment
stat_ihh_class <- R6Class("stat_ihh", inherit = sumstat_class,
  private = list(
    req_segsites = TRUE,
    position = NA,
    population = NULL,
    max_snps = Inf,
    empty_matrix = matrix(0, 0, 3,
                          dimnames = list(NULL, c("IHHa", "IHHd", "IES"))),
    get_snp = function(positions, locus, model) {
      if (is.na(private$position)) return(seq(along = positions))
      pos <- conv_middle_to_trio_pos(private$position, model, locus,
                                     relative_out = FALSE)
      which.min(abs(pos - positions))
    }),
  public = list(
    initialize = function(name, population, position, max_snps) {
      assert_that(is.numeric(population))
      assert_that(length(population) == 1)
      assert_that(is.numeric(max_snps))
      assert_that(length(max_snps) == 1)
      private$population <- population
      private$position <- position
      private$max_snps <- max_snps
      super$initialize(name)
    },
    calculate = function(seg_sites, trees, files, model) {
      assert_that(is.list(seg_sites))
      assert_that(is.model(model))
      pos <- get_snp_positions(seg_sites, model, relative = FALSE)
      ind <- get_population_indiviuals(model, private$population)
      lapply(1:length(seg_sites), function(locus) {
        assert_that(is.matrix(seg_sites[[locus]]))
        if (ncol(seg_sites[[locus]]) == 0) {
          return(private$empty_matrix)
        }
        snps <- private$get_snp(pos[[locus]], locus, model)
        rehh_data <- self$create_rehh_data(seg_sites[[locus]],
                                           pos[[locus]],
                                           ind)
        if (!is.na(private$position)) {
          assert_that(length(snps) == 1)
          ihh <- matrix(0, 1, 3)
          colnames(ihh) <- c("IHHa", "IHHd", "IES")
          ihh[1, 1:2] <- calc_ehh(rehh_data, mrk = snps, plotehh = FALSE)$ihh #nolint
          ihh[1, 3] <- calc_ehhs(rehh_data, mrk = snps, plotehhs = FALSE)$ies #nolint
          return(ihh)
        }
        rehh::scan_hh(rehh_data)[ , -(1:3), drop = FALSE] #nolint
      })
    },
    create_rehh_data = function(seg_sites, pos, ind) {
      assert_that(is.matrix(seg_sites))
      snp_mask <- self$create_snp_mask(seg_sites)
      rehh_data <- new("haplohh")
      rehh_data@haplo <- seg_sites[ind, snp_mask, drop = FALSE] + 1
      rehh_data@position <- pos[snp_mask]
      rehh_data@snp.name <- as.character(seq(along = rehh_data@position))
      rehh_data@chr.name <- 1
      rehh_data@nhap <- length(ind)
      rehh_data@nsnp <- length(rehh_data@position)
      rehh_data
    },
    create_snp_mask = function(seg_sites) {
      n_snps <- ncol(seg_sites)
      if (n_snps < private$max_snps) return(rep(TRUE, n_snps))
      sample.int(n_snps, private$max_snps, replace = FALSE)
    }
  )
)


#' Integrated Extended Haplotype Homozygosity
#'
#' This summary statistic calculates a modified version of the iHH
#' and iES statistics introduced by
#'
#'  Voight et al., A map of recent positive selection in the human genome.
#'  PLoS Biol, 4(3):e72, Mar 2006
#'
#' Coala relies on \code{\link[rehh]{scan_hh}} from package \pkg{rehh} to
#' calculate this statistic. Please refer
#' to their documentation for detailed information on the implementation.
#'
#' @inheritParams sumstat_four_gamete
#' @param position A position relative to the locus extent, e.g. 0.5 for the
#'   middle of the locus. If provided, the statistic will be calculate
#'   for the SNP closest to the given position.
#'   Otherwise, it will be calculated for all SNPs.
#'   The position is relative to the middle locus" extend if trios
#'   are used.
#' @param max_snps The maximal number of SNPs per locus that are used for the
#'   calculation. If a locus has more SNPs than this number, only a
#'   evenly distributed subset of this size will be used for calculating iHS
#'   to increase performance. Set to \code{Inf} to use all SNPs.
#' @return When added to a model, the statistic returns a matrix for each locus.
#'   The columns of the values state the values for integrated EHH for the
#'   ancestral allele (IHHa), integrated EHH for the derived allele (IHHd) and
#'   integrated EHHS (IES), either for all SNPs when no position is given or
#'   for the SNP nearest to the selected position. Each SNP is represented by
#'   a row, sorted by position on the locus.
#' @export
sumstat_ihh <- function(name = "ihh", position = NA, population = 1,
                        max_snps = 1000) {
  stat_ihh_class$new(name, population, position, max_snps)
}
