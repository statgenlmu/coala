#include <fstream>
#include "../inst/include/coala.h"


using namespace Rcpp;

// Reads the 'positions:' line of ms' output and converts it into an
// NumericVector
// [[Rcpp::export]]
NumericVector parse_ms_positions(const std::string line) {
  if (line.size() < 11 || line.substr(0, 11) != "positions: ") {
    Rcout << "Error: Failed to parse positions from ms\' output" << std::endl
          << "Offending line: '" << line << "'" << std::endl;
    stop("Failed to read positions from ms' output");
  }

  std::stringstream stream(line);
  std::vector<double> data;

  // Remove the 'positions: ' at the line's beginning
  stream.ignore(11);

  // Convert the positions into doubles
  std::copy(std::istream_iterator<double>(stream),
  std::istream_iterator<double>(),
  std::back_inserter(data));

  return(wrap(data));
}

// Reads one more files with ms output, and generates a list of segregating
// sites
// [[Rcpp::export]]
List parse_ms_output(const List file_names,
                     const NumericVector sample_size,
                     const int loci_number) {

  std::string line;
  size_t individuals = sum(sample_size);

  List seg_sites(loci_number);
  List trees(loci_number);
  CharacterVector locus_trees;

  int locus = -1;

  for (int i = 0; i < file_names.size(); ++i) {
    CharacterVector file_name = as<CharacterVector>(file_names(i));
    for (int j = 0; j < file_name.size(); ++ j) {

      // Open the file
      std::ifstream output(as<std::string>(file_name(j)).c_str(),
                           std::ifstream::in);

      if (!output.good()) {
        stop(std::string("Cannot open file ") + file_name(0));
      }

      std::getline(output, line);
      // Read it line by line and read the relevant parts
      while (output.good()) {
         // Rcout << "Line: " << line << std::endl;
        if (line == "//") {
          ++locus;
          if (locus >= loci_number) stop("Too many loci in ms output");
          // Rcout << "Locus: " << locus << std::endl;
          std::getline(output, line);
        }

        // Parse Segregating Sites
        else if (line.substr(0, 9) == "segsites:") {
          // Rcout << "Parsing Seg. Sites" << std::endl;
          if (line.substr(0, 11) == "segsites: 0") {
            seg_sites[locus] =
              coala::createSegsites(NumericMatrix(individuals, 0),
                                    NumericVector(0));
          } else {
            std::getline(output, line);

            // Parse Seg.Sites
            NumericVector positions = parse_ms_positions(line);
            NumericMatrix ss(individuals, positions.size());

            for (size_t i = 0; i < individuals; ++i) {
              std::getline(output, line);
              for (int j = 0; j < positions.size(); ++j) {
                ss(i,j) = (line[j] == '1');
              }
            }

            seg_sites[locus] =
              coala::createSegsites(ss, positions, NumericVector(0), false);
          }
          std::getline(output, line);
        }

        // Parse Trees
        else if (line.substr(0, 1) == "[" || line.substr(0, 1) == "(") {
          // Rcout << "Parsing Trees" << std::endl;
          while (line.substr(0, 1) == "[" || line.substr(0, 1) == "(") {
            locus_trees.push_back(line);
            std::getline(output, line);
          }
          trees[locus] = locus_trees;
          locus_trees = CharacterVector();
        }

        else std::getline(output, line);
      }
    }
  }

  if (locus != loci_number - 1) stop("Too few loci in ms output");

  return List::create(_["seg_sites"] = seg_sites,
                      _["trees"] = trees);
}
