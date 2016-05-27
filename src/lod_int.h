// calculate lod support intervals
#ifndef LOD_INT_H
#define LOD_INT_H

#include <Rcpp.h>

// this is the "plain" version ignoring the possibility of multiple LOD peaks
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a pair of indexes (in 0, 1, 2, ..., lod.size()-1) with endpoints of the interval
// followed by all the indexes where LOD == maximum
//
// The R_ version is a wrapper for R
//
Rcpp::IntegerVector R_lod_int_plain(const Rcpp::NumericVector &lod,
                                    const double drop);

std::vector<int> lod_int_plain(const Rcpp::NumericVector& lod,
                               const double drop);

// here we know the peak position and we're looking within a contained subinterval (left, right)
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a pair of indexes (in 0, 1, 2, ..., lod.size()-1) with endpoints of the interval
// followed by all the indexes where LOD == maximum
//
std::vector<int> lod_int_contained(const Rcpp::NumericVector& lod,
                                   const double peakindex, // index in (0,1,2, ..., n-1) where n = lod.size()
                                   const double drop,
                                   const int start,
                                   const int end);

// now the version to deal with a chromosome with multiple peaks
//
// well, first a version with a given peak
// input is
//     lod       : vector of lod scores
//     peakindex : index (in 0,1,...,n-1) of peak
//     drop      : amount to drop for support interval
//     peakdrop  : amount to drop between peaks
// this really only makes sense if peakdrop > drop
// output is just like lod_int_plain
//
std::vector<int> lod_int_peak(const Rcpp::NumericVector& lod,
                              const double peakindex,
                              const double peakdrop,
                              const double drop);

#endif // LOD_INT_H
