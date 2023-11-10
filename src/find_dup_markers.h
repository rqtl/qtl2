// find subsets of markers with identical genotypes
#ifndef FIND_DUP_MARKERS_H
#define FIND_DUP_MARKERS_H

#include <Rcpp.h>

Rcpp::IntegerVector find_dup_markers_notexact(const Rcpp::IntegerMatrix& Geno, // matrix of genotypes, individuals x markers
                                              const Rcpp::IntegerVector& order, // vector indicating order to be considered, most data to least
                                              const Rcpp::IntegerVector markerloc, // integer vector indicating "position"
                                              const bool adjacent_only);   // if true, consider only adjacent markers

#endif // FIND_DUP_MARKERS_H
