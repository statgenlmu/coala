#include <cmath>
#include "../inst/include/coala.h"

using namespace Rcpp;


void maxsplit(const coala::SegSites segsites,
              const int trio_locus,
              const NumericVector individuals,
              const int ploidy,
              int & max_number,
              int & snp_number,
              double & tree) {

  NumericVector trio_locus_vec = coala::getTrioLocus(segsites);
  NumericMatrix ss = coala::getSNPs(segsites);

  std::map<unsigned int, unsigned int> m;
  std::map<unsigned int, std::vector<unsigned int> > balance;
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

    std::vector<unsigned int> sample(individuals.size());
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

	    sample[k] = genotype;
	    key += genotype;
	  }

	  //Associate the key with the corresponding genotype
	  std::map<unsigned int, std::vector<unsigned int> >::iterator b=balance.find(key);
	  if(b==balance.end()) {
	    balance[key]=sample;
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
  unsigned int iden = 0;
  for (std::map<unsigned int,unsigned int>::iterator p=m.begin(); p!=m.end(); ++p) {
	  if(p->second > maxi) {
	    maxi = p->second;
	    iden = p->first;
	  }
  }
  max_number += maxi;

  //Return the genotype associated to the most common split using the key
  std::map<unsigned int, std::vector<unsigned int> >::iterator search = balance.find(iden);
  int count = 0;
  int ones = 0;
  if(search != balance.end()) {
    std::vector<unsigned int> summary = search->second;
    for (std::vector<unsigned int>::iterator p = summary.begin(); p != summary.end(); ++p) {
      ones += *p;
      count++;
    }
  }

  if (count != 0 ) {
    tree = (double)ones/count;
  }else{
    tree = NA_REAL;
  }
}

// [[Rcpp::export]]
NumericMatrix calc_mcmf(const List seg_sites,
                        const NumericVector individuals,
                        const NumericMatrix locus_length,
                        const bool improved = true,
                        const bool has_trios = true,
                        const int ploidy = 1) {

  size_t n_loci = seg_sites.size();

  //Create the matrix that contains the mcmf
  int col_num = 3;
  if (!improved) col_num = 1;

  NumericMatrix mcmf(n_loci, col_num);

  if (improved) {
    mcmf.attr("dimnames") =
      List::create(R_NilValue, CharacterVector::create(
          "mcmf", "bal", "perc_polym")
    );
  }else{
    mcmf.attr("dimnames") =
      List::create(R_NilValue, CharacterVector::create(
          "mcmf")
      );
  }

  coala::SegSites ss;
  int max_split = 0, snp_number = 0, ignore_result = 0;
  double weigth = 0;

  for (size_t locus = 0; locus < n_loci; ++locus) {
    ss = as<coala::SegSites>(seg_sites[locus]);
    NumericVector trio_locus_v;
    trio_locus_v = coala::getTrioLocus(ss);

    if (max(individuals) * ploidy > coala::getSNPs(ss).nrow()) {
      stop("Invalid individuals");
    }

    max_split = 0;
    snp_number = 0;

    weigth = 0;
    if (has_trios) {
      maxsplit(ss, -1, individuals, ploidy, max_split, snp_number, weigth);
      maxsplit(ss, 1, individuals, ploidy, max_split, snp_number, weigth);
      maxsplit(ss, 0, individuals, ploidy, ignore_result, snp_number, weigth);
    } else {
      maxsplit(ss, 0, individuals, ploidy, max_split, snp_number, weigth);
    }

    if (snp_number == 0) {
      mcmf(locus,0) = NA_REAL;
      if(improved) {
        mcmf(locus,1) = NA_REAL;
        mcmf(locus,2) = NA_REAL;
      }
      continue;
    }
    mcmf(locus,0) = (double)max_split / snp_number;
    if(improved) {
      if (weigth < 0.3) {
        weigth = 0;
      }else{
        weigth = 1;
      }
      mcmf(locus, 1) = weigth;

      // Calculate SNPs per basepair
      mcmf(locus, 2) = sum(trio_locus_v == 0) / locus_length(0, 2);
    }
  }

  return(mcmf);
}
