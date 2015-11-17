#include "../inst/include/coala.h"

using namespace Rcpp;


//' Segregating Sites
//'
//' This functions allow the creation and modification of segregating sites
//' objects, which are one of the most basic intermediary statistics that is
//' calculated in coala. Segregating sites are primarily SNP matrix that
//' contains all SNPs for one locus, with some additional information attached.
//' The parts of the S3 class are detailed below.
//'
//' A segregating sites object contains all SNPs for one genetic locus. Each
//' object consists of tree parts: A SNP matrix, a vector of SNP positons and
//' a vector that states which transcript a SNP belong to, if the locus
//' consists of multiple transscripts ('locus trio').
//'
//' \itemize{
//'   \item{In the \strong{SNP} matrix, each row represents a haplotype and each
//'         column represents a SNP. An entry is either \code{1} if the
//'         haplotype carries the derived allele for the SNP, or \code{0} if it
//'         carries the ancestral one.}
//'   \item{In the \strong{positions} vector, each entry gives the relative
//'         position of SNP in the corresponding column of the SNP matrix.}
//'   \item{The \strong{trio_locus} vector contains the trio locus each SNP
//'         belongs to. Entry of \code{-1},\code{0}, \code{1} represent the
//'         left, middle, and right locus, respectively. For normal loci,
//'         this just consists of \code{0}'s}
//' }
//'
//' @param snps The SNP Matrix (see Details).
//' @param positions A numeric vector indicating the relative positions of each
//'   SNP on the locus (see Details).
//' @param trio_locus If the locus consists of a locus trio (see Details).
//' @param check Whether non-segregating sites are remove from the segregating
//'   sites (\code{TRUE}) or not (\code{FALSE}).
//' @export
//' @aliases segsites
//' @author Paul Staab
//'
// [[Rcpp::export]]
coala::SegSites create_segsites(NumericMatrix snps,
                                NumericVector positions,
                                NumericVector trio_locus = NumericVector(0),
                                bool check = true) {
  return coala::createSegsites(snps, positions, trio_locus, check);
}


//' @describeIn create_segsites Returns the SNP matrix from a segregating sites
//'  object.
//'
//' @param segsites The segregating sites object
//' @export
// [[Rcpp::export]]
NumericMatrix get_snps(const coala::SegSites segsites) {
  return coala::getSNPs(segsites);
}


//' @describeIn create_segsites Returns the SNP's positions from a segregating
//'   sites  object.
//' @export
// [[Rcpp::export]]
NumericVector get_positions(const coala::SegSites segsites) {
  return coala::getPositions(segsites);
}


//' @describeIn create_segsites Sets the SNP's positions in a segregating
//'   sites object.
//' @export
// [[Rcpp::export]]
coala::SegSites set_positions(coala::SegSites segsites,
                              const NumericVector positions) {
  return coala::setPositions(segsites, positions);
}


//' @describeIn create_segsites Returns the trio locus positions from a
//'   segregating sites  object.
//' @export
// [[Rcpp::export]]
NumericVector get_trio_locus(const coala::SegSites segsites) {
  return coala::getTrioLocus(segsites);
}


//' @describeIn create_segsites Sets the trio locus in a segregating sites
//'   object.
//' @export
// [[Rcpp::export]]
coala::SegSites set_trio_locus(coala::SegSites segsites,
                               const NumericVector trio_locus) {
  return coala::setTrioLocus(segsites, trio_locus);
}
