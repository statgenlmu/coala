#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;


void addTree(const std::string &tree,
             const size_t length,
             const size_t locus,
             CharacterVector &left,
             CharacterVector &middle,
             CharacterVector &right) {

  if (length == 0) return;
  if (locus % 2 == 1) return;

  std::stringstream ss;
  ss << length;
  std::string prefix = std::string("[") + ss.str() + std::string("]");

  if (locus == 0) {
    left.push_back(prefix + tree);
  } else if (locus == 2) {
    middle.push_back(prefix + tree);
  } else {
    right.push_back(prefix + tree);
  }
}


// [[Rcpp::export]]
List generate_trio_trees(const List trees,
                         const NumericMatrix llm) {

  CharacterVector locus_trees;
  std::string tree;
  size_t digits, pos, locus, len, seg_len, locus_end;

  size_t i = 0, llm_n_row = llm.nrow();
  List result = List(llm_n_row);

  for (size_t llm_row = 0; llm_row < llm_n_row; ++llm_row) {
    CharacterVector left;
    CharacterVector middle;
    CharacterVector right;

    for (size_t llm_row_counter = 0; llm_row_counter < llm(llm_row, 5); ++llm_row_counter) {
      locus_trees = as<CharacterVector>(trees(i++));

      digits = 0;            // Number of digits of the length of the tree
      pos = 0;               // Current position on the locus trio
      locus = 0;             // The locus that we are currently in
      len = 0;               // The length of the current tree
      seg_len = 0;
      locus_end = llm(llm_row, 0); // The end of the current locus

      for (int j = 0; j < locus_trees.size(); ++j) {
        tree = as<std::string>(locus_trees[j]);

        // Get the number of bases for which the tree is valid
        if (tree.substr(0, 1) != "[") {
          len = llm(llm_row, 0) +
                llm(llm_row, 1) +
                llm(llm_row, 2) +
                llm(llm_row, 3) +
                llm(llm_row, 4);
        } else {
          digits = tree.find("]")-1;
          len = std::atoi(tree.substr(1, digits).c_str());
          tree = tree.substr(digits+2, std::string::npos);
        }

        // If the current tree is valid for a sequence that ends behind the
        // end of the locus, we need to do a few things:
        while (pos + len >= locus_end) {
          // First print a tree spanning until the end of the current locus
          seg_len = locus_end - pos;
          addTree(tree, seg_len, locus, left, middle, right);

          // Now the position move towards the end of the current locus.
          len -= seg_len;
          pos += seg_len;

          // And look at the next locus.
          if (locus < 4) {
            ++locus;
            locus_end += llm(llm_row, locus);
          } else {
            // The last tree should end exactly end at the end of the last locus
            if (len != 0) stop("Tree and locus length do not match.");
            pos = 0;

            // Reset the stats and go one to the next locus trio.
            locus = 0;
            locus_end = llm(llm_row, 0);
          }
        }

        // If we are here the tree should end within the current locus.
        // Print the rest and move until to the end of the tree.
        if (len > 0) {
          pos += len;
          addTree(tree, len, locus, left, middle, right);
        }
      }
    }

    if (pos != 0) stop("Error parsing trees");
    result[llm_row] = List::create(_["left"] = left,
                                   _["middle"] = middle,
                                   _["right"] = right);
  }

  return(result);
}
