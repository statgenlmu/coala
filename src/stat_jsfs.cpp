#include "../inst/include/coala.h"

using namespace Rcpp;

//' Calculates the JSFS for two populations
//'
//' @param seg_sites_list List of segregating sites
//' @param pop1 The rows of \code{seg_sites} that correspond to individuals
//'   of the first population.
//' @param pop2 same as \code{pop1}, but for the second population.
//' @export
//' @return The Joint Site Frequency Spectrum, as a matrix.
// [[Rcpp::export]]
NumericMatrix calc_jsfs(const List seg_sites_list,
                        const NumericVector pop1,
                        const NumericVector pop2) {

  size_t n_pop1 = pop1.size();
  size_t n_pop2 = pop2.size();
  size_t nrows_required;
  if (n_pop2 == 0) nrows_required = max(pop1);
  else nrows_required = max(NumericVector::create(max(pop1), max(pop2)));

  NumericMatrix jsfs(n_pop1+1, n_pop2+1);
  size_t idx1, idx2, ncol, nrow;
  coala::SegSites segsites;
  NumericMatrix snps;
  NumericVector trio_locus;

  for (int locus = 0; locus < seg_sites_list.size(); ++locus) {
    segsites = as<coala::SegSites>(seg_sites_list[locus]);
    trio_locus = coala::getTrioLocus(segsites);
    snps = coala::getSNPs(segsites);
    ncol = snps.ncol();
    nrow = snps.nrow();

    if (ncol == 0) continue;
    if (nrow < nrows_required) stop("Seg. Sites has too few rows.");

    for (size_t j = 0; j < ncol; ++j) {
      if (trio_locus(j) != 0) continue; // Only calculate for middle locus

      idx1 = 0; idx2 = 0;

      for (size_t i = 0; i < n_pop1; ++i) idx1 += snps(pop1[i]-1,j);
      for (size_t i = 0; i < n_pop2; ++i) idx2 += snps(pop2[i]-1,j);

      //Rcout << "SNP " << j << ": " << idx1 << "-" << idx2 << std::endl;
      ++jsfs(idx1, idx2);
    }
  }

  return jsfs;
}
