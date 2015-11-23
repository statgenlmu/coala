#include "../inst/include/coala.h"

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calc_nucleotide_div(const ListOf<coala::SegSites> segsites_list,
                                  const NumericVector individuals) {

  size_t n_loci = segsites_list.size();
  NumericVector pi(n_loci);

  NumericMatrix ss;
  size_t n = individuals.size(), m;
  double cnt, c = 2.0 / (n * (n - 1));

  for (size_t locus = 0; locus < n_loci; ++locus) {
    ss = coala::getSNPs(segsites_list[locus]);
    m = ss.ncol();
    cnt = 0;

    for (size_t k = 0; k < m; ++k) {
      for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
          cnt += (ss(individuals[i] - 1, k) != ss(individuals[j] - 1, k));
        }
      }
    }

    pi[locus] = cnt * c;
  }

  return pi;
}
