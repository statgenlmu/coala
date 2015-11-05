#include <RcppArmadilloExtensions/sample.h>
#include "../inst/include/coala.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List unphase_segsites(const List seg_sites_list,
                      const long unsigned int ploidy,
                      const long unsigned int samples_per_ind) {

  if (seg_sites_list.size() == 0) return(List());
  coala::SegSites segsites = as<coala::SegSites>(seg_sites_list[0]);
  NumericMatrix target, source = coala::getSNPs(segsites);
  size_t individuals = source.nrow() / ploidy;
  size_t locus_number = seg_sites_list.size();
  size_t target_rows = individuals * samples_per_ind;
  IntegerVector idxs, possible_idxs = seq_len(ploidy) - 1;
  List ret(locus_number);
  size_t source_offset, target_offset;

  for (size_t locus = 0; locus < locus_number; ++locus) {
    segsites = as<coala::SegSites>(seg_sites_list[locus]);
    source = coala::getSNPs(segsites);
    size_t n_snps = source.ncol();
    target = NumericMatrix(target_rows, n_snps);

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

    ret[locus] = coala::createSegsites(target,
                                       coala::getPositions(segsites),
                                       coala::getTrioLocus(segsites),
                                       ploidy != samples_per_ind);
  }

  return(ret);
}
