#include <RcppArmadillo.h>
#include "../inst/include/coala.h"

using namespace Rcpp;

NumericMatrix read_sequence(CharacterVector output,
                            int &line_nr,
                            const int individuals,
                            const unsigned int locus_length) {

  NumericMatrix seq(individuals, locus_length);
  std::string line, tmp;
  size_t seq_nr, pos;

  for (int i = 0; i < individuals; ++i) {
    ++line_nr;

    // Get the sequence number
    line = output(line_nr);
    tmp = line.substr(0, 10);
    // ms from phyclust adds an "s" to the number
    if (tmp.compare(0, 1, "s") == 0) tmp.erase(0, 1);
    seq_nr = atoi(tmp.c_str()) - 1;

    // Get the sequence
    pos = 0;
    while (pos < locus_length) {
      if (pos == 0) {
        tmp = line.substr(10, std::string::npos);
      } else {
        tmp = output(++line_nr);
      }

      // Translate the positions to numeric values and write them to the matrix
      for (std::string::size_type j = 0; j < tmp.size(); ++j) {
        if (tmp[j] == 'A') seq(seq_nr, pos + j) = 1;
        else if (tmp[j] == 'C') seq(seq_nr, pos + j) = 2;
        else if (tmp[j] == 'G') seq(seq_nr, pos + j) = 3;
        else if (tmp[j] == 'T') seq(seq_nr, pos + j) = 4;
        else {
          /*
          Rcerr << "Error parsing seq-gen sequence " << seq_nr + 1
                << " position " << pos + j + 1 << " of " << locus_length
                << std::endl
                << "Character " << tmp[j] << std::endl
                << "Seq: " << tmp << std::endl;
          */
          stop(std::string("unexpected sequence character: ") + tmp[j]);
        }
      }

      pos += tmp.size();
    }
    if (pos != locus_length) stop("Unexpected locus length.");
  }

  seq.attr("levels") = CharacterVector::create("A", "C", "G", "T");
  return(seq);
}


List conv_seq_to_segsites(NumericMatrix sequence,
                          const int outgroup_size) {

  if (outgroup_size < 1) stop("Outgroup needed to calculate seg. sites");

  int locus_length = sequence.ncol();
  int individuals = sequence.nrow();

  // Determine which positions are SNPs
  std::vector<double> positions;
  int derived_count;
  int outgroup_count;

  for (int j = 0; j < locus_length; ++j) {
    // ignore SNPs with a polymorphic outgroup
    outgroup_count = 0;
    for (int i=individuals-outgroup_size; i<individuals-1; ++i) {
      outgroup_count += (sequence(i, j) != sequence(individuals-1, j));
    }
    if (outgroup_count > 0) continue;

    derived_count = 0;
    for (int i=0; i<individuals-outgroup_size; ++i) {
      derived_count += (sequence(i, j) != sequence(individuals-1, j));
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
        seg_sites(i, derived_count) = (sequence(i, *it) != sequence(individuals-1, *it));
        ++derived_count;
      }
    }

    for (std::vector<double>::iterator it = positions.begin(); it != positions.end(); ++it) {
      *it /= (locus_length - 1);
    }
  }

  return coala::createSegsites(seg_sites, wrap(positions),
                               NumericVector(0), false);
}


// [[Rcpp::export]]
List parse_seqgen_output(CharacterVector output,
                         const int individuals,
                         const int locus_length,
                         const int locus_number,
                         const int outgroup_size,
                         const bool calc_segsites = true) {

  List results(locus_number);
  std::string line;
  int locus = -1;
  NumericMatrix sequence;

  for (int line_nr = 0; line_nr < output.size(); ++line_nr) {
    line = output(line_nr);

    if (line.at(0) == ' ') {
      ++locus;
      if (locus == locus_number) stop("More loci than expected");
      sequence = read_sequence(output, line_nr,
                               individuals, locus_length);
      if (calc_segsites) {
        results(locus) = conv_seq_to_segsites(sequence, outgroup_size);
      } else {
        results(locus) = sequence;
      }

    } else {
      stop(std::string("Unexpect line in seqgen output: ") + line);
    }
  }

  if (locus != locus_number - 1) stop("Less loci than expected");

  return(results);
}
