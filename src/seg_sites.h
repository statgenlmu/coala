#ifndef jaatha_src_seg_sites
#define jaatha_src_seg_sites

namespace Rcpp {
  
inline NumericVector getPositions(NumericMatrix seg_sites) {
  if (!seg_sites.hasAttribute("positions")) stop("SegSites without positions");
  return seg_sites.attr("positions");
}

inline NumericVector getTrioLocus(NumericMatrix seg_sites) {
  if (!seg_sites.hasAttribute("locus")) {
    return(NumericVector(seg_sites.ncol(), 0.0));
  }
  return seg_sites.attr("locus");
}
  
}

#endif