#' @importFrom R6 R6Class
#' @importClassesFrom rehh haplohh
#' @importFrom rehh calc_ehh calc_ehhs ihh2ihs
#' @importFrom methods new checkAtAssignment
stat_ihh_class <- R6Class("stat_ihh", inherit = sumstat_class,
  private = list(
    req_segsites = TRUE,
    position = NA,
    population = NULL,
    max_snps = Inf,
    use_ihs = FALSE,
    empty_matrix = matrix(0, 0, 3,
                          dimnames = list(NULL, c("IHHa", "IHHd", "IES"))),
    get_snp = function(positions, locus, model) {
      if (is.na(private$position)) return(seq(along = positions))
      pos <- conv_middle_to_trio_pos(private$position, model, locus,
                                     relative_out = FALSE)
      which.min(abs(pos - positions))
    }),
  public = list(
    initialize = function(name, population, position, max_snps, calc_ihs) {
      assert_that(is.numeric(population))
      assert_that(length(population) == 1)
      assert_that(is.numeric(max_snps))
      assert_that(length(max_snps) == 1)
      assert_that(is.logical(calc_ihs))
      assert_that(length(calc_ihs) == 1)
      private$population <- population
      private$position <- position
      private$max_snps <- max_snps
      private$use_ihs <- calc_ihs
      super$initialize(name)
    },
    calculate = function(seg_sites, trees, files, model) {
      assert_that(is.list(seg_sites))
      assert_that(is.model(model))
      pos <- get_snp_positions(seg_sites, model, relative = FALSE)
      ind <- get_population_indiviuals(model, private$population)
      lapply(1:length(seg_sites), function(locus) {
        assert_that(is.matrix(seg_sites[[locus]]))
        if (ncol(seg_sites[[locus]]) == 0) return(private$empty_matrix)

        snps <- private$get_snp(pos[[locus]], locus, model)
        rehh_data <- self$create_rehh_data(seg_sites[[locus]],
                                           pos[[locus]],
                                           ind)
        if (rehh_data@nsnp == 0) return(private$empty_matrix)

        ihh <- rehh::scan_hh(rehh_data)

        if (private$use_ihs) {
          if ((rehh_data@nsnp < 50)) freqbin <- 0.90
          else if ((rehh_data@nsnp < 100)) freqbin <- 0.45
          else if ((rehh_data@nsnp < 200)) freqbin <- 0.225
          else if ((rehh_data@nsnp < 400)) freqbin <- 0.1
          else freqbin <- 0.05
          ihs <- suppressWarnings(ihh2ihs(ihh, freqbin))
          ihh <- cbind(ihh, iHS = ihs$res.ihs[ , "iHS"])
        }

        if (!is.na(private$position)) {
          assert_that(length(snps) == 1)
          ihh <- ihh[snps, , drop = FALSE]
        }

        ihh[ , -(1:3), drop = FALSE]
      })
    },
    create_rehh_data = function(seg_sites, pos, ind) {
      assert_that(is.matrix(seg_sites))
      snp_mask <- self$create_snp_mask(seg_sites, ind)
      rehh_data <- new("haplohh")
      rehh_data@haplo <- seg_sites[ind, snp_mask, drop = FALSE] + 1
      rehh_data@position <- pos[snp_mask]
      rehh_data@snp.name <- as.character(seq(along = rehh_data@position))
      rehh_data@chr.name <- 1
      rehh_data@nhap <- length(ind)
      rehh_data@nsnp <- length(rehh_data@position)
      rehh_data
    },
    create_snp_mask = function(seg_sites, ind) {
      polym_in_sample <- apply(seg_sites[ind, , drop = FALSE], 2, function(x) {
        any(0 == x) & any(1 == x)
      })
      n_snps <- sum(polym_in_sample)
      if (n_snps < private$max_snps) return(polym_in_sample)
      sample(which(polym_in_sample), private$max_snps, replace = FALSE)
    }
  )
)


#' Integrated Extended Haplotype Homozygosity
#'
#' This summary statistic calculates a modified version of the iHH,
#' iES and optionally iHS statistics introduced by
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
#'   random subset will be used for calculating iHS to increase performance.
#'   Set to \code{Inf} to use all SNPs.
#' @param calc_ihs If set to \code{TRUE}, additionally standardized iHS is
#'   calculated.
#' @return When added to a model, the statistic returns a matrix for each locus.
#'   The columns of the values contain the values for integrated EHH for the
#'   ancestral allele (IHHa), integrated EHH for the derived allele (IHHd),
#'   integrated EHHS (IES) and iHS (only if \code{calc_ihs = TRUE}),
#'   either for all SNPs when no position is given or
#'   for the SNP nearest to the selected position. Each SNP is represented by
#'   a row, sorted by position on the locus.
#' @export
sumstat_ihh <- function(name = "ihh", position = NA, population = 1,
                        max_snps = 1000, calc_ihs = FALSE) {
  stat_ihh_class$new(name, population, position, max_snps, calc_ihs)
}
