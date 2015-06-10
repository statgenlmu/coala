#' @importFrom R6 R6Class
SumstatIhh <- R6Class('sumstat_ihh', inherit = Sumstat, #nolint
  private = list(
    req_segsites = TRUE,
    position = NA,
    population = NULL,
    empty_matrix = matrix(0, 0, 3,
                          dimnames = list(NULL, c("IHHa", "IHHd", "IES"))),
    get_snp = function(positions, locus, model) {
      if (is.na(private$position)) return(seq(along = positions))
      pos <- conv_middle_to_trio_pos(private$position, model, locus,
                                     relative_out = FALSE)
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
          return(private$empty_matrix)
        }
        snps <- private$get_snp(pos[[locus]], locus, model)
        haplohh <- self$segsites_to_rehh_data(seg_sites[[locus]],
                                              pos[[locus]],
                                              ind)
        if (!is.na(private$position)) {
          assert_that(length(snps) == 1)
          ihh <- matrix(0, 1, 3)
          colnames(ihh) <- c("IHHa", "IHHd", "IES")
          ihh[1, 1:2] <- rehh::calc_ehh(haplohh, mrk = snps, plotehh = FALSE)$ihh
          ihh[1, 3] <- rehh::calc_ehhs(haplohh, mrk = snps, plotehh = FALSE)$ies
          return(ihh)
        }
        rehh::scan_hh(haplohh)[snps , -(1:3), drop = FALSE] #nolint
      })
    },
    segsites_to_rehh_data = function(seg_sites, pos, ind) {
      rehh_data <- new("haplohh")
      rehh_data@haplo <- seg_sites[ind, ] + 1
      rehh_data@position <- pos
      rehh_data@snp.name <- as.character(seq(along = pos))
      rehh_data@chr.name <- 1
      rehh_data@nhap <- length(ind)
      rehh_data@nsnp <- ncol(seg_sites)
      rehh_data
    }
  )
)

#' Integrated Extended Haplotype Homozygosity
#'
#' This summary statistic calculates a modified version of the iHH statistic
#' and iES introduced by
#'
#'  Voight et al., A map of recent positive selection in the human genome.
#'  PLoS Biol, 4(3):e72, Mar 2006
#'
#' Coala relies on \code{\link[rehh]{scan_hh}} from package \pkg{rehh} to
#' calculate this statistic. Please refer
#' to their documentation for detailed information on the concrete
#' implementation.
#' It is required to install the package \pkg{rehh} to use this
#' function.
#'
#' @inheritParams sumstat_four_gamete
#' @param position A position relative to the locus extent, e.g. 0.5 for the
#'   middle of the locus. If provided, the statistic will be calculate
#'   for the SNP closest to the given position.
#'   Otherwise, it will be calculated for all SNPs.
#'   The position is relative to the middle locus' extend if trios
#'   are used.
#' @return When added to a model, the statistic returns a matrix for each locus.
#'   The columns of the values state the values for integrated EHH for the
#'   ancestral allele (IHHa), integrated EHH for the derived allele (IHHd) and
#'   integrated EHHS (IES), either for all SNPs when no position is given or
#'   for the SNP nearest to the selected position. Each SNP is represented by
#'   a row, sorted by position on the locus.
#' @export
sumstat_ihh <- function(name = 'ihh', position=NA, population=1) {
  SumstatIhh$new(name, position, population) #nolint
}
