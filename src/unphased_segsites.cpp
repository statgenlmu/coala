#include <RcppArmadilloExtensions/sample.h>
#include "seg_sites.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List unphase_segsites(const List seg_sites,
                      const long unsigned int ploidy,
                      const long unsigned int samples_per_ind) {

  NumericMatrix target, source = as<NumericMatrix>(seg_sites[0]);
  size_t individuals = source.nrow() / ploidy;
  size_t locus_number = seg_sites.size();
  size_t target_rows = individuals * samples_per_ind;
  IntegerVector idxs, possible_idxs = seq_len(ploidy) - 1;
  List ret(locus_number);
  size_t derived_cnt, new_n_snps, source_offset, target_offset;

  for (size_t locus = 0; locus < locus_number; ++locus) {
    source = as<NumericMatrix>(seg_sites[locus]);
    size_t n_snps = source.ncol();
    target = NumericMatrix(target_rows, n_snps);
    NumericVector positions = getPositions(source);
    NumericVector trio_locus = getTrioLocus(source);

    // Sample pseudo chromosomes
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

    // Remove rows that are no longer polymorphic
    if (ploidy != samples_per_ind) {
      LogicalVector is_snp = no_init(n_snps);
      for (size_t snp = 0; snp < n_snps; ++snp) {
        derived_cnt = sum(target(_, snp));
        is_snp[snp] = !(derived_cnt == 0 || derived_cnt == target_rows);
      }

      new_n_snps = sum(is_snp);
      NumericMatrix out = no_init(target_rows, new_n_snps);
      NumericVector pos = no_init(new_n_snps);
      NumericVector loc = no_init(new_n_snps) ;

      for (size_t i = 0, j = 0; i < n_snps; i++) {
        if (is_snp[i]) {
          out(_,j) = target(_,i);
          pos[j] = positions[i];
          loc[j] = trio_locus[i];
          j = j + 1;
        }
      }

      target = out;
      target.attr("positions") = pos;
      target.attr("locus") = loc;
    } else {
      target.attr("positions") = positions;
      target.attr("locus") = trio_locus;
    }

    ret[locus] = target;
  }

  return(ret);
}
