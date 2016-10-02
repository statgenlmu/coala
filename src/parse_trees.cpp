#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;


void addTree(const std::string &tree,
             const size_t length,
             const size_t trio_locus,
             std::ofstream &left,
             std::ofstream &middle,
             std::ofstream &right) {

  if (length == 0) return;
  if (trio_locus % 2 == 1) return;

  std::stringstream ss;
  ss << length;
  std::string prefix = std::string("[") + ss.str() + std::string("]");

  if (trio_locus == 0) {
    left << prefix + tree << std::endl;
  } else if (trio_locus == 2) {
    middle << prefix + tree << std::endl;
  } else {
    right << prefix + tree << std::endl;
  }
}


// [[Rcpp::export]]
CharacterVector generate_trio_trees(const List trees,
                                    const NumericVector trio_dists,
                                    const CharacterVector file_names) {

  if (trio_dists.size() != 5) stop("trio_dists needs to be of length 5");
  if (file_names.size() != 3) stop("file_names needs to contain 3 files");

  std::ofstream file_left(file_names[0]);
  std::ofstream file_middle(file_names[1]);
  std::ofstream file_right(file_names[2]);

  CharacterVector locus_trees;
  std::string tree;
  size_t digits, pos, locus, len, seg_len, locus_end;
  size_t n_loci = trees.size();

  for (size_t i = 0; i < n_loci; ++i) {
    locus_trees = as<CharacterVector>(trees(i));

    digits = 0;            // Number of digits of the length of the tree
    pos = 0;               // Current position on the locus trio
    locus = 0;             // The locus that we are currently in
    len = 0;               // The length of the current tree
    seg_len = 0;
    locus_end = trio_dists[0]; // The end of the current locus

    for (int j = 0; j < locus_trees.size(); ++j) {
      tree = as<std::string>(locus_trees[j]);

      // Get the number of bases for which the tree is valid
      if (tree.substr(0, 1) != "[") {
        len = sum(trio_dists);
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
        addTree(tree, seg_len, locus, file_left, file_middle, file_right);

        // Now the position move towards the end of the current locus.
        len -= seg_len;
        pos += seg_len;

        // And look at the next locus.
        if (locus < 4) {
          ++locus;
          locus_end += trio_dists[locus];
        } else {
          // The last tree should end exactly end at the end of the last locus
          if (len != 0) stop("Tree and locus length do not match.");
          pos = 0;

          // Reset the stats and go one to the next locus trio.
          locus = 0;
          locus_end = trio_dists[0];
        }
      }

      // If we are here the tree should end within the current locus.
      // Print the rest and move until to the end of the tree.
      if (len > 0) {
        pos += len;
        addTree(tree, len, locus, file_left, file_middle, file_right);
      }
    }

    if (pos != 0) stop("Error parsing trees");
  }

  return(file_names);
}
