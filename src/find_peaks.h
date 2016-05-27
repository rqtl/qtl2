// find peaks in a lod curve for one chromosome
// above threshold, and dropping by at least peakdrop between peaks
#ifndef FIND_PEAKS_H
#define FIND_PEAKS_H

#include <Rcpp.h>

// this version returns a single index for each peak
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a vector of indexes (in 0, 1, 2, ..., lod.size()-1) of peak locations
//
std::vector<int> find_peaks_plain(const Rcpp::NumericVector& lod,
                                  const double threshold,
                                  const double peakdrop);

// like the plain version, but also returning the locations of the valleys in-between
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a vector of two vectors of indexes (in 0, 1, 2, ..., lod.size()-1)
//    - peak locations
//    - valleys between peaks (including 0 and n-1)
//
std::vector< std::vector<int> > find_peaks_valleys(const Rcpp::NumericVector& lod,
                                                   const double threshold,
                                                   const double peakdrop);

// this version deals with ties in the LOD scores
// ...it returns all indexes which jointly achieve the maximum LOD score
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a list of vectors of indexes (in 0, 1, 2, ..., lod.size()-1) of peak locations
//
std::vector< std::vector<int> > find_peaks(const Rcpp::NumericVector& lod,
                                           const double threshold,
                                           const double peakdrop);


// this version returns both peaks and lod intervals
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a list of vectors of indexes (in 0, 1, 2, ..., lod.size()-1)
//     first two values are the left and right endpoints of the interval
//     remaining values are the indexes with the maximum LOD score
//
std::vector< std::vector<int> > find_peaks_and_lodint(const Rcpp::NumericVector& lod,
                                                      const double threshold,
                                                      const double peakdrop,
                                                      const double drop);

#endif // FIND_PEAKS_H
