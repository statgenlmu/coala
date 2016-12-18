#include <cmath>
#include "../inst/include/coala.h"

using namespace Rcpp;
using namespace std;

void maxsplit(const coala::SegSites segsites,
              const int trio_locus,
              const NumericVector individuals,
              const int ploidy,
              int & max_number,
              int & snp_number,
              double & weigth) {

  NumericVector trio_locus_vec = coala::getTrioLocus(segsites);
  NumericMatrix ss = coala::getSNPs(segsites);

  std::map<unsigned int, unsigned int> m;
  std::map<unsigned int, std::vector<unsigned int> > balance;
  unsigned int key, ind_nr, genotype;
  bool switch_snp;
  bool switch_snp_decided;
  double half_ploidy = 0.5 * ploidy;

  for (int snp = 0; snp < ss.ncol(); ++snp) { //Loop thorugh the columns
    if (trio_locus_vec[snp] != trio_locus) continue;

    // For each SNP, create a binary representation of the SNP...
    switch_snp = false;
    switch_snp_decided = false;
    key = 0;

    std::vector<unsigned int> sample(individuals.size());
	  for (int k = 0; k < individuals.size(); ++k) { //Loop through the rows
	    key *= ploidy + 1;

      ind_nr = (individuals[k] - 1) * ploidy;
	    genotype = 0;
	    for (int chr = 0; chr < ploidy; ++chr) genotype += ss(ind_nr + chr, snp); //ss(2row, 1column)

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

    // Ignore snps that are not segregating (in at least one of) the individuals
    if (key == 0) continue;

	  // and increase the corresponding counter
	  ++snp_number; //This is K in the document, be careful because the other have been left out
	  std::map<unsigned int,unsigned int>::iterator p=m.find(key); //Look for the key
	  if(p==m.end()) { //If the key doesn't exist, the counter will reach the end of p
	    m[key]=1; //Then initialize the counter for that key
	  } else {
	    ++(p->second); //Otherwise, if the key already exists, increase the counter one time
	  }
  }

  // Return the maximal counter
  unsigned int maxi = 0;
  unsigned int iden = 0; //Get the key of the maximal counter
  for (std::map<unsigned int,unsigned int>::iterator p=m.begin(); p!=m.end(); ++p) {
	  if(p->second > maxi) { //Get the value in m whose counter achieved the largest number
	    maxi = p->second;
	    iden = p->first; //and also get the key for the next step
	  }
  }
  max_number += maxi;

  //Return the genotype associated to the most common split using the key(iden)
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
    weigth = (double)ones/count;
  }else{
    weigth = NA_REAL;
  }
}

// [[Rcpp::export]]
NumericMatrix calc_mcmf(const List seg_sites,
                        const NumericVector individuals,
                        const bool has_trios = true,
                        const bool expand_mcmf = false,
                        const int type_expand = 1,
                        const int ploidy = 1,
                        const NumericMatrix locus_length = NumericMatrix(0)) {

  size_t n_loci = seg_sites.size();

  //Create the matrix that contains the mcmf
  int col_num;
  if (type_expand == 1) {
    col_num = 1;
  } else if (type_expand == 2) {
    col_num = 2;
  } else {
    col_num = 3;
  }

  NumericMatrix mcmf(n_loci, col_num);

  if (type_expand == 1 ) {
    mcmf.attr("dimnames") =
      List::create(R_NilValue, CharacterVector::create(
          "mcmf")
      );
  }else if( type_expand == 2 ) {
    mcmf.attr("dimnames") =
      List::create(R_NilValue, CharacterVector::create(
          "mcmf", "bal")
      );
  }else{
    mcmf.attr("dimnames") =
      List::create(R_NilValue, CharacterVector::create(
          "mcmf", "bal", "perc_polym")
      );
  }

  coala::SegSites ss;
  int max_split = 0, snp_number = 0, ignore_result = 0;
  double weigth = 0;

  for (size_t locus = 0; locus < n_loci; ++locus) { //Loop through the loci
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
      if(type_expand == 2 || type_expand == 3) {
        mcmf(locus,1) = NA_REAL;
      }
      if(type_expand == 3) {
        mcmf(locus,2) = NA_REAL;
      }
      continue;
    }

    mcmf(locus,0) = (double)max_split / snp_number;

    if (type_expand == 2 || type_expand == 3) {
      mcmf(locus, 1) = weigth;
      if (type_expand == 3) {
        // Calculate SNPs per basepair
        mcmf(locus, 2) = sum(trio_locus_v == 0) / locus_length(0, 2);
      }
    }
  }

  return(mcmf);
}
