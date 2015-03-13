#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;


// [[Rcpp::export]]
List parse_trees(const List file_names,
                 const int loci_number,
                 const bool separate_loci = true) {

  int file_number = file_names.size();
  int locus = -1;
  int file_nr = 0;

  List trees;
  if (separate_loci) trees = List(loci_number);
  else trees = List(file_number);

  CharacterVector locus_trees;
  CharacterVector file_name;
  std::string line;

  for (int i = 0; i < file_names.size(); ++i) {
    file_name = as<CharacterVector>(file_names(i));
    if (file_name.size() != 1) stop("Expecting one file per locus");

    // Open the file
    std::ifstream input(as<std::string>(file_name(0)).c_str(),
                        std::ifstream::in);
    if (!input.is_open()) {
      stop(std::string("Cannot open file ") + file_name(0));
    }

    while (getline(input, line)) {
      if (line.substr(0, 2) == "//") {
        if (separate_loci) {
          if (locus >= 0) trees(locus) = locus_trees;
          ++locus;
          locus_trees = CharacterVector();
        }
      }
      if (line.substr(0, 1) == "[" || line.substr(0, 1) == "(") {
        locus_trees.push_back(line);
      }
    }
    if (!separate_loci) {
      trees(file_nr) = locus_trees;
      ++file_nr;
      locus_trees = CharacterVector();
    }
  }

  if (separate_loci) trees(locus) = locus_trees;

  return(trees);
}


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

  if (trees.size() != llm.nrow()) stop("Locus number mismatch.");

  CharacterVector locus_trees;
  std::string tree;
  size_t digits, pos, locus, len, seg_len, locus_end;
  List result = List(trees.size());


  for (size_t i = 0; i < trees.size(); ++i) {
    locus_trees = as<CharacterVector>(trees[i]);

    digits = 0;            // Number of digits of the length of the tree
    pos = 0;               // Current position on the locus trio
    locus = 0;             // The locus that we are currently in
    len = 0;               // The length of the current tree
    seg_len = 0;
    locus_end = llm(i, 0); // The end of the current locus

    CharacterVector left;
    CharacterVector middle;
    CharacterVector right;

    for (size_t j = 0; j < locus_trees.size(); ++j) {
      tree = as<std::string>(locus_trees[j]);

      // Get the number of bases for which the tree is valid
      digits = tree.find("]")-1;
      len = std::atoi(tree.substr(1, digits).c_str());
      tree = tree.substr(digits+2, std::string::npos);

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
          locus_end += llm(i, locus);
          //if (locus == 2) current_tree = &trio_trees[1];
          //if (locus == 4) current_tree = &trio_trees[2];
        } else {
          // The last tree should end exactly end at the end of the last locus
          if (len != 0) stop("Tree and locus length do not match.");
          // Reset the stats and go one to the next locus trio.
          locus = 0;
          locus_end = llm(i, 0);
          pos = 0;
          //current_tree = &trio_trees[0];
        }
      }

      // If we are here the tree should end within the current locus.
      // Print the rest and move until to the end of the tree.
      if (len > 0) {
        pos += len;
        addTree(tree, len, locus, left, middle, right);
      }
    }

    if (pos != 0) stop("Error parsing trees");
    result[i] = List::create(_["left"] = left,
                             _["middle"] = middle,
                             _["right"] = right);
  }

  return(result);
}
