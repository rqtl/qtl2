// pick subset of well-spaced markers
//
// More preciesly, find subset of markers that maximizes sum(weights)
// subject to the condition that no two adjacent markers are within
// distance d.

#ifndef PICK_MARKER_SUBSET_H
#define PICK_MARKER_SUBSET_H

#include <Rcpp.h>

Rcpp::IntegerVector pick_marker_subset(const Rcpp::NumericVector& pos,       // positions of markers
                                       const double min_d,                   // minimum position between markers
                                       const Rcpp::NumericVector& weights);  // weights on the markers

#endif // PICK_MARKER_SUBSET_H
