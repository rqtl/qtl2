// Utilities to interpolate between genetic and physical maps
#ifndef INTERPOLATE_MAPS_H
#define INTERPOLATE_MAPS_H

#include <Rcpp.h>

// find interval in map that contains pos
// [-1 if to left, map.size()-1 if to right]
// map should be sorted
int find_interval(const double pos, const Rcpp::NumericVector& map);

// for positions relative to oldmap, interpolate to get positions relative to newmap
Rcpp::NumericVector interpolate_map(const Rcpp::NumericVector& oldpos,
                                    const Rcpp::NumericVector& oldmap,
                                    const Rcpp::NumericVector& newmap);

#endif // INTERPOLATE_MAPS_H
