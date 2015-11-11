#include "../inst/include/coala.h"

using namespace Rcpp;


//' Segregating Sites
//'
//' This functions allow the creation and modification of segregating sites
//' objects, which are on of the most basic statistics that is calculated in
//' coala. Segregating Sites are primarily a matrix, where each row repesents
//' a haplotype and each column represents a SNPs. A given entry is either 1 if
//' the haplotype carries the derived allele for the SNP, or 0 if it carries
//' the ancestral one.
//'
//' @param snps The SNP Matrix
//' @param positions A numeric vector indicating the relative positions of each
//'   SNP on the locus.
//' @param trio_locus If the locus consists of a locus trio, then this
//'   contains the trio locus each SNP belongs to.
//' @param check Whether non-segregating sites are remove from the segregating
//'   sites (\code{TRUE}) or not (\code{FALSE}).
//' @export
//' @aliases segsites
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
