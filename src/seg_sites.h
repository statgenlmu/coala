#ifndef coala_src_seg_sites
#define coala_src_seg_sites

namespace coala {

inline Rcpp::NumericMatrix
removeFixedPos(const Rcpp::NumericMatrix seg_sites);

inline Rcpp::NumericMatrix
createSegsites(Rcpp::NumericMatrix snps,
               Rcpp::NumericVector positions = Rcpp::NumericVector(0),
               Rcpp::NumericVector trio_locus = Rcpp::NumericVector(0),
               bool check = true) {

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

  // Check that there are no fixed positions in the seg sites.
  if (check) snps = coala::removeFixedPos(snps);

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

inline Rcpp::NumericMatrix
removeFixedPos(const Rcpp::NumericMatrix seg_sites) {
  size_t n_snps = seg_sites.ncol(),
         n_ind = seg_sites.nrow();
  size_t derived_cnt, new_n_snps;

  Rcpp::LogicalVector is_snp = Rcpp::no_init(n_snps);
  for (size_t snp = 0; snp < n_snps; ++snp) {
    derived_cnt = Rcpp::sum(seg_sites(Rcpp::_, snp));
    is_snp[snp] = !(derived_cnt == 0 || derived_cnt == n_ind);
  }

  new_n_snps = Rcpp::sum(is_snp);
  if (new_n_snps == n_snps) return(seg_sites);

  Rcpp::NumericMatrix out = Rcpp::no_init(n_ind, new_n_snps);
  Rcpp::NumericVector pos = Rcpp::no_init(new_n_snps);
  Rcpp::NumericVector loc = Rcpp::no_init(new_n_snps) ;
  Rcpp::NumericVector positions = getPositions(seg_sites);
  Rcpp::NumericVector trio_locus = getTrioLocus(seg_sites);

  for (size_t i = 0, j = 0; i < n_snps; i++) {
    if (is_snp[i]) {
      out(Rcpp::_,j) = seg_sites(Rcpp::_,i);
      pos[j] = positions[i];
      loc[j] = trio_locus[i];
      j = j + 1;
    }
  }

  return createSegsites(out, pos, loc, false);
}

}

#endif
