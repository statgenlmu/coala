#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;

// Function that reads a file containing output from ms or msms ('in_file')
// and copies the trees in newick format into a separate file ('out_file').
// This new file can be used as input for seq-gen. Only Unix system, this could
// be done with a few grep's, but this function should give us some platform 
// independence.
// [[Rcpp::export]]
CharacterVector parseTrees(const std::string in_file, 
                           const NumericVector trio_opts,
                           Function tempfile) {
    
  std::string line; 
  bool trio = false;
  if (trio_opts(0) != 0 ||  trio_opts(4) != 0) trio = true;
  
  // Prepare a file for the first locus' tree
  std::vector<std::string> out_files;
  out_files.push_back(as<std::string>(tempfile("trees")));

  // Open input file and first output file
  std::ifstream input(in_file.c_str(), std::ifstream::in);
  if (!input.is_open()) 
    stop(std::string("Failed to open simulation results in ") + in_file);
  
  std::ofstream output_0(out_files[0].c_str(), std::ofstream::out);
  if (!output_0.is_open()) stop("Failed to create a file for writing trees.");
  std::ofstream* output = &output_0;

  // Filter trees
  if (!trio) {
    // Fast parser when not using loci-trios.
    while (getline(input, line)) {
      if (line.substr(0, 1) == "[")  *output << line << "\n";
    }
  } else {
    // Open remaining output files
    out_files.push_back(as<std::string>(tempfile("trees")));
    out_files.push_back(as<std::string>(tempfile("trees")));
    
    std::ofstream output_1(out_files[1].c_str(), std::ofstream::out);
    if (!output_1.is_open()) stop("Failed to create a file for writing trees.");
    std::ofstream output_2(out_files[2].c_str(), std::ofstream::out);
    if (!output_2.is_open()) stop("Failed to create a file for writing trees.");
    
    // When using loci-trios, we need to remove the inter-locus regions.
    size_t digits = 0,  // Number of digits of the length of the tree
    pos = 0,     // Current position on the locus trio
    locus = 0,   // The locus that we are currently in
    len = 0,     // The length of the current tree
    seg_len = 0, 
    locus_end = trio_opts[0]; // The end of the current locus
    
    while (getline(input, line)) {
      // Read in the next tree
      if (line.substr(0, 1) != "[") continue;
      digits = line.find("]")-1;
      len = std::atoi(line.substr(1, digits).c_str());
      
      // If the current tree is valid for a sequence that ends behind the 
      // end of the locus, we need to do a few things:
      while (pos + len >= locus_end) {
        // First print a tree spanning until the end of the current locus 
        // if neccessary.
        seg_len = locus_end - pos;
        //Rprintf("locus %i (end %i) - pos %i - len %i\n", locus, locus_end, 
        //        pos, len);
        
        if (locus % 2 == 0) {
          *output << "[" << seg_len << "]" 
                  << line.substr(digits+2, std::string::npos) << "\n";
        }
        
        // Now the position move towards the end of the current locus.
        len -= seg_len;
        pos += seg_len;
        
        // And look at the next locus.
        if (locus < 4) {
          ++locus;
          locus_end += trio_opts[locus];
          if (locus == 2) output = &output_1;
          if (locus == 4) output = &output_2;
        } else {
          // The last tree should end exactly end at the end of the last locus
          if (len != 0) stop("Tree and locus length do not match.");
          // Reset the stats and go one to the next locus trio.
          locus = 0;
          locus_end = trio_opts[0];
          pos = 0;
          output = &output_0;
        }
      }
      
      // If we are here the tree should end within the current locus.
      // Print the rest and move until to the end of the tree.
      if (len > 0) {
        pos += len;
        if (locus % 2 == 0) {
          *output << "[" << len << "]" 
                  << line.substr(digits+2, std::string::npos) << "\n";
        }
      }
    }
    
    if (pos != 0) stop("Error parsing trees");
    output_1.close();
    output_2.close();
  }
  
  input.close();
  output_0.close();

  return(wrap(out_files));
}