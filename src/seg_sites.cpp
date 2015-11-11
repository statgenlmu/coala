#include "../inst/include/coala.h"

using namespace Rcpp;

// [[Rcpp::export]]
coala::SegSites create_segsites(NumericMatrix snps,
                                NumericVector positions = NumericVector(0),
                                NumericVector trio_locus = NumericVector(0),
                                bool check = true) {
  return coala::createSegsites(snps, positions, trio_locus, check);
}

// [[Rcpp::export]]
NumericMatrix get_snps(const coala::SegSites seg_sites) {
  return coala::getSNPs(seg_sites);
}


// [[Rcpp::export]]
NumericVector get_positions(const coala::SegSites seg_sites) {
  return coala::getPositions(seg_sites);
}

// [[Rcpp::export]]
coala::SegSites set_positions(coala::SegSites seg_sites,
                              const NumericVector positions) {
  return coala::setPositions(seg_sites, positions);
}

// [[Rcpp::export]]
NumericVector get_trio_locus(const coala::SegSites seg_sites) {
  return coala::getTrioLocus(seg_sites);
}

// [[Rcpp::export]]
coala::SegSites set_trio_locus(coala::SegSites seg_sites,
                               const NumericVector trio_locus) {
  return coala::setTrioLocus(seg_sites, trio_locus);
}
