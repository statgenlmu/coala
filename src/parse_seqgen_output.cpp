#include <Rcpp.h>
#include <fstream>

#include "seg_sites.h"

using namespace Rcpp;

NumericMatrix parseSeqgenSegSites(std::ifstream &output,
                                  const int individuals,
                                  const int locus_length,
                                  const int outgroup_size) {
  
  std::string tmp;
  size_t seq_nr;
  
  // First read the complete locus and save it in `sequence`.
  std::vector<std::vector<char> > sequence(individuals);
  for (int i=0; i<individuals; ++i) {
    // Read sequence number
    if (!output.good()) 
    stop("Unexpeced end of seq-gen file.");
    output >> tmp;
    // ms from phyclust adds an "s" to the seqName. Remove it if there.
    if (tmp.compare(0, 1, "s") == 0) tmp.erase(0,1);
    seq_nr = atoi(tmp.c_str());
    
    // Read sequence
    output >> tmp;
    const char* cstr=tmp.c_str();
    sequence.at(seq_nr-1).assign(cstr, cstr+tmp.length());
  }
  
  // Determine which positions are SNPs
  std::vector<double> positions;
  int derived_count;
  int outgroup_count;
  
  for (int j=0; j<locus_length; ++j) {
    // ignore SNPs with a polymorphic outgroup
    outgroup_count = 0;
    for (int i=individuals-outgroup_size; i<individuals-1; ++i) {
      outgroup_count += (sequence.at(i).at(j) != sequence.at(individuals-1).at(j));
    }
    if (outgroup_count > 0) continue;
    
    derived_count = 0;
    for (int i=0; i<individuals-outgroup_size; ++i) {
      derived_count += (sequence.at(i).at(j) != sequence.at(individuals-1).at(j));
    }
    if (derived_count > 0 && derived_count < (individuals - outgroup_size)) {
      positions.push_back(j);
    }
  }
  
  NumericMatrix seg_sites(individuals-outgroup_size, positions.size());
  
  if (positions.size() > 0) {
    for (int i=0; i<individuals-outgroup_size; ++i) {
      derived_count = 0;
      for (std::vector<double>::iterator it = positions.begin(); it != positions.end(); ++it) {
        seg_sites(i, derived_count) = (sequence[i][*it] != sequence[individuals-1][*it]);
        ++derived_count;
      }
    }
    

    for (std::vector<double>::iterator it = positions.begin(); it != positions.end(); ++it) {
        *it /= (locus_length - 1);
    }
  }
  
  seg_sites.attr("positions") = wrap(positions); 
  return seg_sites;
}

NumericMatrix cbindPos(NumericMatrix ss_l, 
                       NumericMatrix ss_m, 
                       NumericMatrix ss_r) {
                         
  size_t snps = ss_l.ncol() + ss_m.ncol() + ss_r.ncol();
  int sample_size = ss_m.nrow();
  NumericMatrix ss(sample_size, snps);
  NumericVector positions(snps);
  NumericVector positions_tmp;
  NumericVector locus(snps);
  
  positions_tmp = getPositions(ss_l);
  for (int col = 0; col < ss_l.ncol(); ++col) {
    for (int row = 0; row < sample_size; ++row) {
      ss(row, col) = ss_l(row, col);
    }
    locus(col) = -1;
    positions(col) = positions_tmp(col);
  }
  
  int offset = ss_l.ncol();
  positions_tmp = getPositions(ss_m);
  for (int col = 0; col < ss_m.ncol(); ++col) {
    for (int row = 0; row < sample_size; ++row) {
      ss(row, col + offset) = ss_m(row, col);
    }
    locus(col + offset) = 0;
    positions(col + offset) = positions_tmp(col);
  }
          
  offset = offset + ss_m.ncol();
  positions_tmp = getPositions(ss_r);
  for (int col = 0; col < ss_r.ncol(); ++col) {
    for (int row = 0; row < sample_size; ++row) {
      ss(row, col + offset) = ss_r(row, col);
    }
    locus(col + offset) = 1;
    positions(col + offset) = positions_tmp(col);
  }
  
  ss.attr("locus") = locus;
  ss.attr("positions") = positions;
  return ss;
}

// [[Rcpp::export]]
List parseSeqgenOutput(const List file_names, 
                       const int sample_size,
                       const NumericMatrix sequence_length,
                       const int loci_number,
                       const int outgroup_size = 1) {

  std::string line_l, line_m, line_r;

  List seg_sites(loci_number);
  int locus = -1;
  NumericVector positions(0);
  
  for (int i = 0; i < file_names.size(); ++i) {
    CharacterVector file_name = as<CharacterVector>(file_names(i));

    // Open the file
    if (file_name.size() == 1) {
      std::ifstream output_m(as<std::string>(file_name(0)).c_str(), std::ifstream::in);
      if (!output_m.is_open()) stop(std::string("Cannot open file ") + file_name(0));

      while ( output_m.good() ) {
        std::getline(output_m, line_m);
        if (line_m == "") continue;
        if (line_m.substr(0, 1) == " ") {
          ++locus;
          seg_sites[locus] = parseSeqgenSegSites(output_m, sample_size, 
                                                 sequence_length(i, 2),
                                                 outgroup_size);
        
        } else {
          stop(std::string("Unexpected line in seq-gen output: '") + line_m + "'");
        }
      }
      
    } else if (file_name.size() == 3) {
      std::ifstream output_l(as<std::string>(file_name(0)).c_str(), std::ifstream::in);
      std::ifstream output_m(as<std::string>(file_name(1)).c_str(), std::ifstream::in);
      std::ifstream output_r(as<std::string>(file_name(2)).c_str(), std::ifstream::in);
      if (!output_l.is_open()) stop(std::string("Cannot open file ") + file_name(0));
      if (!output_m.is_open()) stop(std::string("Cannot open file ") + file_name(1));
      if (!output_r.is_open()) stop(std::string("Cannot open file ") + file_name(2));
      
      // We already know the information in the first line
      while ( output_m.good() && output_l.good() && output_r.good()) {
        std::getline(output_m, line_m);
        std::getline(output_l, line_l);
        std::getline(output_r, line_r);

        if (line_m == "") continue;
        if (line_m.substr(0, 1) == " ") {
          if (line_l.substr(0, 1) != " ") stop("seq-gen outputs not in sync");
          if (line_r.substr(0, 1) != " ") stop("seq-gen outputs not in sync");
          ++locus;
          
          
          NumericMatrix ss_l = parseSeqgenSegSites(output_l, sample_size, 
                                                   sequence_length(i, 0), outgroup_size);
          NumericMatrix ss_m = parseSeqgenSegSites(output_m, sample_size, 
                                                   sequence_length(i, 2), outgroup_size);
          NumericMatrix ss_r = parseSeqgenSegSites(output_r, sample_size, 
                                                   sequence_length(i, 4), outgroup_size);
                                                   
          seg_sites[locus] = cbindPos(ss_l, ss_m, ss_r);
        } else {
          stop("Unexpected line in seq-gen output.");
        }        
      } 

    } else {
      stop("Unexpected number of seq-gen simulation files");
    }
  }
  return(seg_sites);
}


