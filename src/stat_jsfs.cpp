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
  size_t idx1, idx2, ncol;
  NumericMatrix snps;
  NumericVector trio_locus;

  for (int locus = 0; locus < segsites_list.size(); ++locus) {
    trio_locus = coala::getTrioLocus(segsites_list[locus]);
    snps = coala::getSNPs(segsites_list[locus]);
    ncol = snps.ncol();

    if (ncol == 0) continue;

    for (size_t j = 0; j < ncol; ++j) {
      if (trio_locus(j) != 0) continue; // Only calculate for middle locus

      idx1 = 0; idx2 = 0;

      for (size_t i = 0; i < n_inds[0]; ++i) idx1 += snps(ind_per_pop[0][i]-1,j);
      for (size_t i = 0; i < n_inds[1]; ++i) idx2 += snps(ind_per_pop[1][i]-1,j);

      //cout << "SNP " << j << ": " << idx1 << "-" << idx2 << std::endl;
      ++jsfs(idx1 + idx2 * (n_inds[0] + 1));
    }
  }

  if (n_pops == 2) {
    jsfs.attr("class") = "matrix";
    jsfs.attr("dim") = n_inds + 1;
  }

  return jsfs;
}
