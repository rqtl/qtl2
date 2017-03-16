// identify IBD segments in founder genotypes
#ifndef FOUNDER_IBD_SEG_H
#define FOUNDER_IBD_SEG_H

#include <Rcpp.h>

// find_IBD_segments
// For a pair of individuals on a single chromosome:
//   calculate LOD score for each interval for evidence of IBD vs not
//   find set of non-overlapping intervals with the largest LOD scores
//
// Input:
//   g1  Genotypes of individual 1 (values 1/3, no NAs)
//   g2  Genotypes of individual 2 (values 1/3, no NAs)
//   p   Frequency of genotype 1 at each marker
//   error_prob  Probability of error or mutation at a marker
//   (g1, g2, p all have common length, n = number of markers
//
// Output:
//   matrix with n rows and the following 6 columns
//    - left endpoint (just the values 1...n)
//    - right endpoint (as values in 1...n)
//    - LOD score
//    - number of markers
//    - number of mismatches
//    - 0/1 indicating whether to retain (0 = some overlapping interval has larger LOD score)
//
Rcpp::NumericMatrix find_IBD_segments(const Rcpp::IntegerVector& g1,
                                      const Rcpp::IntegerVector& g2,
                                      const Rcpp::NumericVector& p,
                                      const double error_prob);

#endif // FOUNDER_IBD_SEG_H
