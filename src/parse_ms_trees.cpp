#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;

// [[Rcpp::export]]
List parse_ms_trees(const List files,
                    const int loci_number) {

  int locus = -1;
  List trees;
  trees = List(loci_number);

  CharacterVector locus_trees;
  CharacterVector files_vec;
  std::string line;

  for (int i = 0; i < files.size(); ++i) {
    files_vec = as<CharacterVector>(files(i));
    for (int j = 0; j < files_vec.size(); j++) {
      // Open the file
      std::ifstream input(as<std::string>(files_vec(0)).c_str(),
                          std::ifstream::in);
      if (!input.is_open()) {
        stop(std::string("Cannot open file ") + files_vec(0));
      }

      while (getline(input, line)) {
        // A new locus starts
        if (line.substr(0, 2) == "//") {
          if (locus >= 0) trees(locus) = locus_trees;
          ++locus;
          locus_trees = CharacterVector();
        }

        // Save tree lines
        else if (line.substr(0, 1) == "[" || line.substr(0, 1) == "(") {
          locus_trees.push_back(line);
        }
      }
    }
  }

  // Save last locus tree
  trees(locus) = locus_trees;

  return(trees);
}
