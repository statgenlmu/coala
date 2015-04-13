#include <RcppArmadilloExtensions/sample.h>
#include "seg_sites.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List unphase_segsites(const List seg_sites,
                      const int ploidy,
                      const int samples_per_ind) {

  NumericMatrix target, source = as<NumericMatrix>(seg_sites[0]);
  size_t individuals = source.nrow() / ploidy;
  size_t locus_number = seg_sites.size();
  size_t target_rows = individuals * samples_per_ind;
  IntegerVector idxs, possible_idxs = seq_len(ploidy) - 1;
  List ret(locus_number);

  size_t source_offset, target_offset;

  for (size_t locus = 0; locus < locus_number; ++locus) {
    source = as<NumericMatrix>(seg_sites[locus]);
    size_t n_snps = source.ncol();
    target = NumericMatrix(target_rows, n_snps);

    for (size_t ind = 0; ind < individuals; ++ind) {
      source_offset = ind * ploidy;
      target_offset = ind * samples_per_ind;

      for (size_t snp = 0; snp < n_snps; ++snp) {
        idxs = RcppArmadillo::sample(possible_idxs, samples_per_ind, false);
        for (size_t i = 0; i < samples_per_ind; ++i) {
          target(target_offset + i, snp) = source(source_offset + idxs[i], snp);
        }
      }

    }

    target.attr("positions") = source.attr("positions");
    target.attr("locus") = source.attr("locus");
    ret[locus] = target;
  }

  return(ret);
}