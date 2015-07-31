#include <cmath>
#include <Rcpp.h>
#include "seg_sites.h"
using namespace Rcpp;


void maxsplit(const NumericMatrix ss,
              const int trio_locus,
              const NumericVector individuals,
              int & max_number,
              int & snp_number) {

  NumericVector trio_locus_vec = getTrioLocus(ss);

  std::map<unsigned int, unsigned int> m;
  unsigned int key;
  int max_key_value = std::pow(2, individuals.size());

  for (int snp = 0; snp < ss.ncol(); ++snp) {
    if (trio_locus_vec[snp] != trio_locus) continue;

    // For each SNP, create a binary representation of the SNP...
    key = 0;
	  for(int k=1; k < individuals.size(); ++k) {
	    key *= 2;
	    key += (ss(individuals[k]-1, snp) != ss(individuals[0]-1, snp));
	  }

    // Ignore snps that are not segregating in the given individuals
    if (key == 0 || key == max_key_value) continue;

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
                        const bool has_trios = true) {

  size_t n_loci = seg_sites.size();
  NumericVector mcmf(n_loci);

  NumericMatrix ss;
  int max_split = 0, snp_number = 0, ignore_result = 0;

  for (size_t locus = 0; locus < n_loci; ++locus) {
    ss = as<NumericMatrix>(seg_sites[locus]);
    if (max(individuals) > ss.nrow()) stop("Invalid individuals");

    max_split = 0;
    snp_number = 0;

    if (has_trios) {
      maxsplit(ss, -1, individuals, max_split, snp_number);
      maxsplit(ss, 1, individuals, max_split, snp_number);
      maxsplit(ss, 0, individuals, ignore_result, snp_number);
    } else {
      maxsplit(ss, 0, individuals, max_split, snp_number);
    }

    if (snp_number == 0) {
      mcmf[locus] = NA_REAL;
      continue;
    }
    mcmf[locus] = (double)(max_split) / snp_number;
  }

  return(mcmf);
}
