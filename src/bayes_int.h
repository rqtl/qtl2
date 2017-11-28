// calculate Bayes credible intervals
#ifndef BAYES_INT_H
#define BAYES_INT_H

#include <Rcpp.h>

// this is the "plain" version ignoring the possibility of multiple LOD peaks
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a pair of indexes (in 0, 1, 2, ..., lod.size()-1) with endpoints of the interval
// followed by all the indexes where LOD == maximum
//
// The R_ version is a wrapper for R
//
Rcpp::IntegerVector R_bayes_int_plain(const Rcpp::NumericVector& lod,
                                      const Rcpp::NumericVector& pos,
                                      const double prob);

std::vector<int> bayes_int_plain(const Rcpp::NumericVector& lod,
                                 const Rcpp::NumericVector& pos,
                                 const double prob);

// here we know the peak position and we're looking within a contained subinterval (left, right)
//
// input is a vector of LOD scores ordered by position along a chromosome plus
// plus the actual locations
//
// output is a pair of indexes (in 0, 1, 2, ..., lod.size()-1) with endpoints of the interval
// followed by all the indexes where LOD == maximum
//
std::vector<int> bayes_int_contained(const Rcpp::NumericVector& lod,
                                     const Rcpp::NumericVector& pos,
                                     const double peakindex, // index in (0,1,2, ..., n-1) where n = lod.size()
                                     const double prob,
                                     const int start,
                                     const int end);

#endif // BAYES_INT_H
