#include "../inst/include/coala.h"

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
                         const size_t snp,
                         const size_t ploidy) {

  size_t mut_count = 0, ind_offset;
  for (size_t k = 0; k < n_ind; ++k) {
    ind_offset = (individuals[k] - 1) * ploidy;
    for (size_t l = 0; l < ploidy; ++l) {
      mut_count += seg_sites(ind_offset + l, snp);
    }
  }

  return (mut_count <= ploidy || mut_count >= (n_ind - 1) * ploidy);
}


// [[Rcpp::export]]
NumericMatrix calc_four_gamete_stat(const ListOf<coala::SegSites> seg_sites_list,
                                    const IntegerVector individuals,
                                    const NumericMatrix locus_length,
                                    const unsigned int ploidy = 1) {

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
  coala::SegSites seg_sites;
  NumericMatrix snps_matrix;
  NumericVector positions, trio_locus;
  std::vector<size_t> snps;
  snps.reserve(1000);
  size_t n_snps, idx_i, idx_j;
  double sample_prob;

  for (size_t locus = 0; locus < loci_number; ++locus) {
    // Reset variables
    std::fill(total_count.begin(), total_count.end(), 0);

    // Get the locus
    seg_sites = seg_sites_list[locus];
    snps_matrix = coala::getSNPs(seg_sites);
    positions = coala::getPositions(seg_sites);
    trio_locus = coala::getTrioLocus(seg_sites);

    // Filter SNPs which are non-polymorpic or singletons in the population
    n_snps = snps_matrix.ncol();
    for (size_t i = 0; i < n_snps; ++i) {
      if (!is_singleton(snps_matrix, individuals, n_ind, i, ploidy)) {
        snps.push_back(i);
      }
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
        size_t n[2], ind_offset;
        for (size_t k = 0; k < n_ind; ++k) {
          ind_offset = (individuals[k] - 1) * ploidy;
          n[0] = 0;
          n[1] = 0;
          for (size_t l = 0; l < ploidy; ++l) {
            n[0] += snps_matrix(ind_offset + l, idx_i);
            n[1] += snps_matrix(ind_offset + l, idx_j);
          }
          if ((n[0] == 0 || n[0] == ploidy) || (n[1] == 0 || n[1] == ploidy)) {
            for (size_t l = 0; l < ploidy; ++l) {
              combinations[2 * snps_matrix(ind_offset + l, idx_i) +
                            snps_matrix(ind_offset + l, idx_j)] = true;
            }
          }
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
