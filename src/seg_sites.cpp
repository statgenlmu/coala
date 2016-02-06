#include "../inst/include/coala.h"

using namespace Rcpp;


//' @describeIn segsites Creates segregating sites
//'
//' @param snps The SNP Matrix (see Details).
//' @param positions A numeric vector indicating the relative positions of each
//'   SNP on the locus (see Details).
//' @param trio_locus If the locus consists of a locus trio (see Details).
//' @param check Whether non-segregating sites are removed from the segregating
//'   sites (\code{TRUE}) or not (\code{FALSE}).
//' @export
//'
// [[Rcpp::export]]
coala::SegSites create_segsites(NumericMatrix snps,
                                NumericVector positions,
                                NumericVector trio_locus = NumericVector(0),
                                bool check = true) {
  return coala::createSegsites(snps, positions, trio_locus, check);
}


//' @describeIn segsites Returns the SNP matrix from a segregating sites
//'  object.
//'
//' @param segsites The segregating sites object
//' @export
// [[Rcpp::export]]
NumericMatrix get_snps(const coala::SegSites segsites) {
  return coala::getSNPs(segsites);
}


//' @describeIn segsites Returns the SNP's positions from a segregating
//'   sites  object.
//' @export
// [[Rcpp::export]]
NumericVector get_positions(const coala::SegSites segsites) {
  return coala::getPositions(segsites);
}


//' @describeIn segsites Sets the SNP's positions in a segregating
//'   sites object.
//' @export
// [[Rcpp::export]]
coala::SegSites set_positions(coala::SegSites segsites,
                              const NumericVector positions) {
  return coala::setPositions(segsites, positions);
}


//' @describeIn segsites Returns the trio locus positions from a
//'   segregating sites  object.
//' @export
// [[Rcpp::export]]
NumericVector get_trio_locus(const coala::SegSites segsites) {
  return coala::getTrioLocus(segsites);
}


//' @describeIn segsites Sets the trio locus in a segregating sites
//'   object.
//' @export
// [[Rcpp::export]]
coala::SegSites set_trio_locus(coala::SegSites segsites,
                               const NumericVector trio_locus) {
  return coala::setTrioLocus(segsites, trio_locus);
}
