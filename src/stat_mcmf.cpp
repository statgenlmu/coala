#include <cmath>
#include "../inst/include/coala.h"

using namespace Rcpp;


void maxsplit(const coala::SegSites segsites,
              const int trio_locus,
              const NumericVector individuals,
              const int ploidy,
              int & max_number,
              int & snp_number) {

  NumericVector trio_locus_vec = coala::getTrioLocus(segsites);
  NumericMatrix ss = coala::getSNPs(segsites);

  std::map<unsigned int, unsigned int> m;
  unsigned int key, ind_nr, genotype;
  bool switch_snp;
  bool switch_snp_decided;
  double half_ploidy = 0.5 * ploidy;

  for (int snp = 0; snp < ss.ncol(); ++snp) {
    if (trio_locus_vec[snp] != trio_locus) continue;

    // For each SNP, create a binary representation of the SNP...
    switch_snp = false;
    switch_snp_decided = false;
    key = 0;
	  for (int k = 0; k < individuals.size(); ++k) {
	    key *= ploidy + 1;

      ind_nr = (individuals[k] - 1) * ploidy;
	    genotype = 0;
	    for (int chr = 0; chr < ploidy; ++chr) genotype += ss(ind_nr + chr, snp);

	    if (!switch_snp_decided) {
	      if (genotype > half_ploidy) {
	        switch_snp = true;
	        switch_snp_decided = true;
	      } else if (genotype < half_ploidy) {
	        switch_snp_decided = true;
	      }
	    }
	    if (switch_snp) genotype = ploidy - genotype;

	    key += genotype;
	  }

    // Ignore snps that are not segregating in the given individuals
    if (key == 0) continue;

	  // and increase the corresponding counter
	  ++snp_number;
	  std::map<unsigned int,unsigned int>::iterator p=m.find(key);
	  if(p==m.end()) {
	    m[key]=1;
	  } else {
	    ++(p->second);
	  }
  }

  // Return the maximal counter
  unsigned int maxi = 0;
  for (std::map<unsigned int,unsigned int>::iterator p=m.begin(); p!=m.end(); ++p) {
	  if(p->second > maxi) maxi = p->second;
  }
  max_number += maxi;
}



// [[Rcpp::export]]
NumericVector calc_mcmf(const List seg_sites,
                        const NumericVector individuals,
                        const bool has_trios = true,
                        const int ploidy = 1) {

  size_t n_loci = seg_sites.size();
  NumericVector mcmf(n_loci);

  coala::SegSites ss;
  int max_split = 0, snp_number = 0, ignore_result = 0;

  for (size_t locus = 0; locus < n_loci; ++locus) {
    ss = as<coala::SegSites>(seg_sites[locus]);
    if (max(individuals) * ploidy > coala::getSNPs(ss).nrow()) {
      stop("Invalid individuals");
    }

    max_split = 0;
    snp_number = 0;

    if (has_trios) {
      maxsplit(ss, -1, individuals, ploidy, max_split, snp_number);
      maxsplit(ss, 1, individuals, ploidy, max_split, snp_number);
      maxsplit(ss, 0, individuals, ploidy, ignore_result, snp_number);
    } else {
      maxsplit(ss, 0, individuals, ploidy, max_split, snp_number);
    }

    if (snp_number == 0) {
      mcmf[locus] = NA_REAL;
      continue;
    }
    mcmf[locus] = (double)max_split / snp_number;
  }

  return(mcmf);
}

