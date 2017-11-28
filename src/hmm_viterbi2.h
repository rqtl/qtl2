// Viterbi algoirthm to find arg max Pr(g | O)
//   where g = sequence of true genotypes and O = observed marker genotypes
// (this version assumes constant is_female and cross_info and pre-calcs the step and emit matrices)
#ifndef HMM_VITERBI2_H
#define HMM_VITERBI2_H

#include <Rcpp.h>

// find most probable sequence of genotypes
Rcpp::IntegerMatrix viterbi2(const Rcpp::String& crosstype,
                             const Rcpp::IntegerMatrix& genotypes, // columns are individuals, rows are markers
                             const Rcpp::IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                             const bool is_X_chr,
                             const bool is_female, // same for all individuals
                             const Rcpp::IntegerVector& cross_info, // same for all individuals
                             const Rcpp::NumericVector& rec_frac,   // length nrow(genotypes)-1
                             const Rcpp::IntegerVector& marker_index, // length nrow(genotypes)
                             const double error_prob);

#endif // HMM_VITERBI2_H
