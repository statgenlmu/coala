#' @importFrom R6 R6Class
SumstatNsl <- R6Class('sumstat_nsl', inherit = SumstatIhh, #nolint
  public = list(
    segsites_to_snp_map = function(seg_sites, pos) {
      map <- data.frame(name = seq(along = pos),
                        chr = 1,
                        pos = seq(along = pos),
                        anc = 0,
                        der = 1)
      file <- tempfile('snp_map')
      write.table(map, file, row.names = FALSE, col.names = FALSE)
      file
    },
    calculate = function(seg_sites, files, model) {
      stat <- super$calculate(seg_sites, files, model)
      lapply(stat, function(x) {
        colnames(x) <- NULL
        x[ , 3]
      })
    }
  )
)


#' Number  of  Segregating  Sites  by  Length
#'
#' This (mis)uses \code{\link[rehh]{scan_hh}} from package \pkg{rehh} to
#' calculate the nSL statistic from
#'
#' Ferrer-Admetlla et al., On Detecting Incomplete Soft or Hard Selective
#' Sweeps Using Haplotype Structure. Mol Biol Evol (2014) 31 (5): 1275-1291.
#' doi:10.1093/molbev/msu077
#'
#' It uses the package \pkg{rehh} to calculate iHS with distances between
#' SNPs measured in SNPs rather than in physical or genetic distance, which
#' should result in nSL.
#' It is required to install the package \pkg{rehh} to use this function.
#'
#' @inheritParams sumstat_ihh
#' @return When added to a model, the statistic returns a vector for each locus,
#'   which lists the values of nSL either for all SNPs when no position is
#'   given or for the SNP nearest to the selected position.
#'   SNPs are sorted by their positions on the locus.
#' @export
sumstat_nSL <- function(name = 'ihh', position=NA, population=1) { #nolint
  SumstatNsl$new(name, position, population) #nolint
}
