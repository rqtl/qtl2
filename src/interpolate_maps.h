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

// apply find_interval() to each of a vector of positions
Rcpp::IntegerVector find_intervals(const Rcpp::NumericVector& pos,
                                   const Rcpp::NumericVector& map);

// For each position, having already figured out the interval it's in,
// figure out whether it is exactly at the left location or not
Rcpp::LogicalVector is_pos_on_map(const Rcpp::NumericVector& pos,
                                  const Rcpp::NumericVector& map,
                                  const Rcpp::IntegerVector& interval, double tol);

#endif // INTERPOLATE_MAPS_H
