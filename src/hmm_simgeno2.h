// simulate genotypes given observed marker data
// (this version assumes constant is_female and cross_info and pre-calcs the step and emit matrices)

#ifndef HMM_SIMGENO2_H
#define HMM_SIMGENO2_H

#include <Rcpp.h>

// simulate genotypes given observed marker data
Rcpp::IntegerVector sim_geno2(const Rcpp::String& crosstype,
                              const Rcpp::IntegerMatrix& genotypes, // columns are individuals, rows are markers
                              const Rcpp::IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                              const bool is_X_chr,
                              const bool is_female, // same for all individuals
                              const Rcpp::IntegerVector& cross_info, // same for all individuals
                              const Rcpp::NumericVector& rec_frac,   // length nrow(genotypes)-1
                              const Rcpp::IntegerVector& marker_index, // length nrow(genotypes)
                              const double error_prob,
                              const int n_draws); // number of imputations

#endif // HMM_SIMGENO2_H
