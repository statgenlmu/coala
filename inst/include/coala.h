#ifndef coala_inst_include_coala
#define coala_inst_include_coala

#include <RcppArmadillo.h>


/*
 * Documentation for the functions here is included in the file
 * `src/seg_sites.cc`, because it seems that Rcpp currently can't Rcpp::export
 * inlined functions directly.
 */


namespace coala {

typedef Rcpp::List SegSites;

inline SegSites removeFixedPos(const SegSites seg_sites);

inline SegSites
createSegsites(const Rcpp::NumericMatrix snps,
               const Rcpp::NumericVector positions = Rcpp::NumericVector(0),
               Rcpp::NumericVector trio_locus = Rcpp::NumericVector(0),
               const bool check = true) {

  // Check positions
  if (positions.size() != snps.ncol()) {
    Rcpp::stop("Number of positions differs from the number of SNPS");
  }

  // Check trio locus
  if (trio_locus.size() == 0) {
    trio_locus = Rcpp::rep(0, positions.size());
  }
  if (trio_locus.size() != snps.ncol()) {
    Rcpp::stop("Length of trio_locus differs from the number of SNPS");
  }

  // Change the class
  SegSites segsites = Rcpp::List::create(Rcpp::_["snps"] = snps,
                                         Rcpp::_["position"] = positions,
                                         Rcpp::_["trio_locus"] = trio_locus);

  segsites.attr("class") = "segsites";

  // Check that there are no fixed positions in the seg sites.
  if (check) segsites = coala::removeFixedPos(segsites);

  return(segsites);
}


inline Rcpp::NumericMatrix getSNPs(const SegSites seg_sites) {
  return seg_sites["snps"];
}


inline Rcpp::NumericVector getPositions(const SegSites seg_sites) {
  return seg_sites["position"];
}


inline SegSites setPositions(SegSites seg_sites,
                             const Rcpp::NumericVector positions) {
  if (positions.size() != getSNPs(seg_sites).ncol()) {
    Rcpp::stop("Number of positions differs from the number of SNPS");
  }
  seg_sites["position"] = positions;
  return seg_sites;
}


inline Rcpp::NumericVector getTrioLocus(const SegSites seg_sites) {
  return seg_sites["trio_locus"];
}


inline SegSites setTrioLocus(SegSites seg_sites,
                             const Rcpp::NumericVector trio_locus) {
  if (trio_locus.size() != getSNPs(seg_sites).ncol()) {
    Rcpp::stop("Length of trio_locus differs from the number of SNPS");
  }
  seg_sites["trio_locus"] = trio_locus;
  return seg_sites;
}


inline SegSites removeFixedPos(const SegSites seg_sites) {
  Rcpp::NumericMatrix snps = seg_sites["snps"];

  size_t n_snps = snps.ncol(),
         n_ind = snps.nrow();
  size_t derived_cnt, new_n_snps;

  Rcpp::LogicalVector is_snp = Rcpp::no_init(n_snps);
  for (size_t snp = 0; snp < n_snps; ++snp) {
    derived_cnt = Rcpp::sum(snps(Rcpp::_, snp));
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
      out(Rcpp::_,j) = snps(Rcpp::_,i);
      pos[j] = positions[i];
      loc[j] = trio_locus[i];
      j = j + 1;
    }
  }

  return createSegsites(out, pos, loc, false);
}

}

#endif
