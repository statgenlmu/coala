#include <Rcpp.h>
#include "seg_sites.h"

using namespace Rcpp;

inline int getType(size_t snp_1, size_t snp_2,
                   const NumericVector positions,
                   const NumericVector trio_locus) {

  if (trio_locus[snp_1] == 0 && trio_locus[snp_2] == 0) {
    return(std::abs(positions[snp_1] - positions[snp_2]) > 0.1);
  }

  if (trio_locus(snp_1) != trio_locus(snp_2)) {
    return 3;
  }

  return 2;
}

inline bool is_singleton(const NumericMatrix seg_sites,
                         const IntegerVector individuals,
                         const size_t n_ind,
                         const size_t snp) {

  size_t mut_count = 0;
  for (size_t i = 0; i < n_ind; ++i) {
    mut_count += seg_sites(individuals[i]-1, snp);
  }

  //return(mut_count <= 1);
  return(mut_count <= 1 || mut_count >= n_ind-1);
}


// [[Rcpp::export]]
NumericMatrix calc_four_gamete_stat(const List seg_sites_list,
                                    const IntegerVector individuals,
                                    const NumericMatrix locus_length) {

  size_t loci_number = seg_sites_list.size();
  size_t n_ind = individuals.size();
  size_t locus_group = 0;
  size_t locus_in_group = 0;
  size_t type; // 0 = near, 1 = far, 2 = on outer locus, 3 = between loci;

  NumericMatrix violations(loci_number, 6);
  violations.attr("dimnames") =
    List::create(R_NilValue, CharacterVector::create(
        "mid_near", "mid_far", "outer", "between", "mid", "perc_polym")
    );

  std::vector<size_t> total_count(4);
  std::vector<bool> combinations(4);
  NumericMatrix seg_sites;
  NumericVector positions, trio_locus;
  std::vector<size_t> snps;
  snps.reserve(1000);
  size_t n_snps, idx_i, idx_j;
  double sample_prob;

  for (size_t locus = 0; locus < loci_number; ++locus) {
    // Reset variables
    std::fill(total_count.begin(), total_count.end(), 0);

    // Get the locus
    seg_sites = as<NumericMatrix>(seg_sites_list[locus]);
    positions = getPositions(seg_sites);
    trio_locus = getTrioLocus(seg_sites);

    // Filter SNPs which are non-polymorpic or singletons in the population
    n_snps = seg_sites.ncol();
    for (size_t i = 0; i < n_snps; ++i) {
      if (!is_singleton(seg_sites, individuals, n_ind, i)) snps.push_back(i);
    }
    n_snps = snps.size();

    // If we have more than 1250 good pairs, sample approximately 1000
    sample_prob = 1000.0 / (0.5 * n_snps * (n_snps - 1));

    // Look at all pairs of SNPs
    for (size_t i = 1; i < n_snps; ++i) {
      for (size_t j = 0; j < i; ++j) {
        if (sample_prob < 0.75 && runif(1)[0] > sample_prob) continue;
        idx_i = snps[i];
        idx_j = snps[j];

        // Ignore SNP pairs between outer loci
        if (std::abs(trio_locus[idx_i] - trio_locus[idx_j]) == 2) continue;

        // Get the type of the SNP pair
        type = getType(idx_i, idx_j, positions, trio_locus);

        // Reset combination counter
        std::fill(combinations.begin(), combinations.end(), false);

        // Count combinations
        for (size_t k = 0; k < n_ind; ++k) {
          combinations[2*seg_sites(individuals[k]-1, idx_i) +
                         seg_sites(individuals[k]-1, idx_j)] = true;
        }

        // If we have all combinations
        if (combinations[0] && combinations[1] &&
            combinations[2] && combinations[3]) {
          ++violations(locus, type);
        }

        ++total_count[type];
      }
    }

    // Calculate %violations for all SNP pairs in middle locus
    violations(locus, 4) = (violations(locus, 0) + violations(locus, 1)) /
                           (total_count[0] + total_count[1]);

    // Calculate %violations for near and far SNP pairs in middle locus,
    // and for pairs between loci.
    for (int i = 0; i < 4; ++i) violations(locus, i) /= total_count[i];

    // Calculate SNPs per basepair for middle locus
    violations(locus, 5) = sum(trio_locus == 0) / locus_length(locus_group, 2);

    // Clean Up
    snps.clear();
    ++locus_in_group;
    if (locus_in_group == locus_length(locus_group, 5)) {
      locus_in_group = 0;
      ++locus_group;
    }
  }

  return violations;
}
