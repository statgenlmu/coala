#ifndef coala_src_seg_sites
#define coala_src_seg_sites

namespace coala {

inline Rcpp::NumericMatrix
createSegsites(Rcpp::NumericMatrix snps,
               Rcpp::NumericVector positions = Rcpp::NumericVector(0),
               Rcpp::NumericVector trio_locus = Rcpp::NumericVector(0)) {

  // Set positions
  if (positions.size() == 0 && snps.hasAttribute("positions")) {
    positions = snps.attr("positions");
  }
  if (positions.size() != snps.ncol()) {
    Rcpp::stop("Number of positions differs from the number of SNPS");
  }
  snps.attr("positions") = positions;

  // Set trio locus
  if (trio_locus.size() == 0 && snps.hasAttribute("trio_locus")) {
    trio_locus = snps.attr("trio_locus");
  }
  if (trio_locus.size() > 0) snps.attr("trio_locus") = trio_locus;

  // Change the class
  snps.attr("class") =
    Rcpp::CharacterVector::create("segsites", "matrix", "array",
                                  "structure", "vector");

  return(snps);
}

inline Rcpp::NumericVector
getPositions(const Rcpp::NumericMatrix seg_sites) {
  if (!seg_sites.hasAttribute("positions")) {
    Rcpp::stop("SegSites without positions");
  }
  return seg_sites.attr("positions");
}

inline Rcpp::NumericVector
getTrioLocus(const Rcpp::NumericMatrix seg_sites) {
  if (!seg_sites.hasAttribute("trio_locus")) {
    return(Rcpp::NumericVector(seg_sites.ncol(), 0.0));
  }
  return seg_sites.attr("trio_locus");
}

}

#endif
