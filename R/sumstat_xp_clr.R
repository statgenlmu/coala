#' @importFrom R6 R6Class
stat_xp_clr_class <- R6Class("stat_xp_clr", inherit = sumstat_class,
 private = list(
   req_segsites = TRUE,
   pop_focal = NULL,
   pop_reference = NULL,
   empty_matrix = data.frame(CHR = numeric(),
                             POSITION = numeric(),
                             FREQ_a = numeric(),
                             IHHa = numeric(),
                             HHd = numeric(),
                             IES = numeric())
 ),
 public = list(
   initialize = function(name, pop_focal, pop_reference, transformation) {
     assert_that(is.numeric(pop_focal) && length(pop_focal) == 1)
     private$pop_focal <- pop_focal

     assert_that(is.numeric(pop_reference) && length(pop_reference) == 1)
     private$pop_reference <- pop_reference

     super$initialize(name, transformation)
   },
   calculate = function(seg_sites, trees, files, model, sim_tasks = NULL) {
     # Create input files
     geno_file_focal <- public$create_geno_file(seg_sites, model, private$pop_focal)
     geno_file_reference <- public$create_geno_file(seg_sites, model, private$pop_focal)

     # Clean-Up
     unlink(c(geno_file_focal, geno_file_reference))
   },
   create_geno_file = function(seg_sites, model, pop) {
     assert_that(is.list(seg_sites))
     assert_that(is.model(model))

     pop_individuals <- get_population_individuals(model, pop)
     geno_file <- tempfile(paste("xp_clr_pop", pop, "geno", sep = "_"))

     for (locus_seg_sites in seg_sites) {
       assert_that(is_segsites(locus_seg_sites))
       geno_data <- t(locus_seg_sites$snps[pop_individuals, , drop = FALSE])
       write.table(geno_data, file = geno_file, append = TRUE, col.names = FALSE, row.names = FALSE)
     }

     geno_file
   }
 )
)


#' Summary Statistic: Integrated Extended Haplotype Homozygosity
#'
#' This summary statistic calculates a number of values based on
#' extended haplotype homozygosity (EHH), including iHH, iES
#' and optionally iHS.
#' Coala relies on \code{\link[rehh]{scan_hh}} from package \pkg{rehh} to
#' calculate this statistic. Please refer
#' to their documentation for detailed information on the implementation.
#' Please cite the corresponding publication (see below) if you use the
#' statistic for a publication.
#'
#' @section References:
#' \itemize{
#'  \item{Mathieu Gautier and Renaud Vitalis, rehh: an R package to detect
#'       footprints of selection in genome-wide SNP data from
#'       haplotype structure. Bioinformatics (2012) 28 (8): 1176-1177
#'       first published online March 7, 2012 doi:10.1093/bioinformatics/bts115}
#'  \item{Voight et al., A map of recent positive selection in the human
#'                        genome. PLoS Biol, 4(3):e72, Mar 2006.}
#'  }
#'
#' @inheritParams sumstat_four_gamete
#' @param max_snps The maximal number of SNPs per locus that are used for the
#'   calculation. If a locus has more SNPs, only a random subset of them will
#'   be used to increase performance. Set to \code{Inf} to use all SNPs.
#' @param calc_ihs If set to \code{TRUE}, additionally standardized iHS is
#'   calculated.
#' @return If \code{calc_ihs = FALSE}, a data.frame with values for
#'   iHH and iES is returned. Otherwise, a list of two data frames are
#'   returned, one for IHH and IES values and the other one for IHS values.
#'
#'   In all `data.frames` rows are SNPs and the columns present the following
#'   values for each SNP:
#'   \itemize{
#'    \item{CHR: The SNP's locus}
#'    \item{Positions: The SNP's absolute position on its locus}
#'    \item{FREQ_a: The SNP's absolute position on its locus}
#'    \item{IHHa: integrated EHH for the ancestral allele}
#'    \item{IHHd: integrated EHH for the derived allele}
#'    \item{IES: integrated EHHS}
#'    \item{iHS: iHS, normalized over all loci.}
#'   }
#' @export
#' @template summary_statistics
#' @examples
#'   model <- coal_model(20, 1, 1000) +
#'     feat_mutation(1000) +
#'     sumstat_ihh()
#' \dontrun{
#'     stat <- simulate(model)
#'     print(stat$ihh)}
#' @author Paul Staab
sumstat_xp_clr <- function(name = "xp_clr", pop_focal, pop_reference,
                           transformation = identity) {
  stat_xp_clr_class$new(name, pop_focal, pop_reference, transformation)
}
