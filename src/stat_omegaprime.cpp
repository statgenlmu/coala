#include <cmath>
#include <Rcpp.h>
#include "seg_sites.h"
using namespace Rcpp;


int maxsplit(const NumericMatrix ss,
             const int trio_locus,
             const NumericVector individuals) {

  NumericVector trio_locus_vec = getTrioLocus(ss);

  std::map<unsigned int,unsigned int> m;
  unsigned int key;

  for (int snp = 0; snp < ss.ncol(); ++snp) {
    if (trio_locus_vec[snp] != trio_locus) continue;

    // For each SNP, create a binary representation of the SNP...
    key=0;
	  for(int k=1; k < individuals.size(); ++k) {
	    key *= 2;
	    key += (ss(individuals[k]-1, snp) != ss(individuals[0]-1, snp));
	  }

    // Ignore snps that are not segregating in individuals
    if (key == 0) continue;

	  // and increase the corresponding counter
	  std::map<unsigned int,unsigned int>::iterator p=m.find(key);
	  if(p==m.end()) {
	    m[key]=1;
	  } else {
	    ++(p->second);
	  }
  }

  // Return the maximal counter
  unsigned int maxi=0;
  for (std::map<unsigned int,unsigned int>::iterator p=m.begin(); p!=m.end(); ++p) {
	  if(p->second > maxi) maxi = p->second;
  }
  return maxi;
}


// [[Rcpp::export]]
NumericVector calc_omegaprime(List seg_sites, NumericVector individuals) {
  size_t n_loci = seg_sites.size();
  NumericVector omega_prime(n_loci);

  NumericMatrix ss;

  for (size_t locus = 0; locus < n_loci; ++locus) {
    ss = as<NumericMatrix>(seg_sites[locus]);
    if (max(individuals) > ss.nrow()) stop("Invalid individuals");
    if (ss.ncol() == 0) omega_prime[locus] = NA_REAL;
    else {
          omega_prime[locus] = (double)(maxsplit(ss, -1, individuals) +
                                  maxsplit(ss, 1, individuals)) / ss.ncol();
    }
  }

  return(omega_prime);
}
