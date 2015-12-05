// reduce markers to a more evenly-spaced set
#ifndef REDUCE_MARKERS_H
#define REDUCE_MARKERS_H

#include <Rcpp.h>

// pos = vector of marker positions
// Seek subset of markers such that no two markers are within min_dist
// and sum(weights) is as large as possible
//
// return value is integer vector of marker indices,
// in {1, 2, ..., length(pos)}
Rcpp::IntegerVector reduce_markers(const Rcpp::NumericVector& pos,
                                   const Rcpp::NumericVector& weights,
                                   const double min_dist);

#endif // REDUCE_MARKERS_H
