#include <Rcpp.h>
#include "seg_sites.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix create_segsites(NumericMatrix snps,
                              NumericVector positions = NumericVector(0),
                              NumericVector trio_locus = NumericVector(0)) {
  return coala::createSegsites(snps, positions, trio_locus);
}


// [[Rcpp::export]]
NumericVector get_positions(const NumericMatrix seg_sites) {
  return coala::getPositions(seg_sites);
}


// [[Rcpp::export]]
NumericVector get_trio_locus(const NumericMatrix seg_sites) {
  return coala::getTrioLocus(seg_sites);
}
