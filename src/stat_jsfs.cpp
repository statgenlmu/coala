#include "../inst/include/coala.h"
#include <cmath>

using namespace Rcpp;

//' Calculates the JSFS
//'
//' @param segsites_list List of segregating sites
//' @param ind_per_pop A list of integer vector, where each entry gives the
//'   index of the haploids that belong the the corresponding population.
//' @export
//' @return The Joint Site Frequency Spectrum, as a matrix.
// [[Rcpp::export]]
NumericVector calc_jsfs(const ListOf<coala::SegSites> segsites_list,
                        const ListOf<IntegerVector> ind_per_pop) {

  size_t n_pops = ind_per_pop.size();

  NumericVector n_inds = no_init(n_pops);
  size_t n_entries = 1;

  for (size_t i = 0; i < n_pops; ++i) {
    n_inds(i) = ind_per_pop[i].size();
    n_entries *= n_inds[i] + 1;
  }

  NumericVector jsfs(n_entries, 0);
  size_t ncol, k;
  NumericMatrix snps;
  NumericVector trio_locus;

  for (int locus = 0; locus < segsites_list.size(); ++locus) {
    trio_locus = coala::getTrioLocus(segsites_list[locus]);
    snps = coala::getSNPs(segsites_list[locus]);
    ncol = snps.ncol();

    if (ncol == 0) continue;

    for (size_t j = 0; j < ncol; ++j) {
      // Only include SNPs on the middle locus
      if (trio_locus[j] != 0) continue;

      // idx will be a multiindex, where idx[i] is a number
      // of derived alleles in population i.
      std::vector<size_t> idx(n_pops, 0);

      for(size_t n=0; n < n_pops; ++n) {
        for(size_t i = 0; i < n_inds[n]; ++i) {
          idx[n] += snps(ind_per_pop[n][i]-1, j);
        }
      }

      // this will be the index in the output jsfs
      k = idx[n_pops-1];

      for(int m = n_pops-2; m >= 0; --m) {
        k *= (n_inds[m]+1);
        k += idx[m];
      }

      ++jsfs[k];
    }
  }

  if (n_pops == 2) {
    jsfs.attr("class") = "matrix";
    jsfs.attr("dim") = n_inds + 1;
  } else if (n_pops > 2) {
    jsfs.attr("class") = "array";
    jsfs.attr("dim") = n_inds + 1;
  }

  return jsfs;
}
