// simulate genotypes given observed marker data
#ifndef HMM_SIMGENO_H
#define HMM_SIMGENO_H

#include <Rcpp.h>

// simulate genotypes given observed marker data
Rcpp::IntegerVector sim_geno(const Rcpp::String& crosstype,
                             const Rcpp::IntegerMatrix& genotypes, // columns are individuals, rows are markers
                             const Rcpp::IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                             const bool is_X_chr,
                             const Rcpp::LogicalVector& is_female, // length n_ind
                             const Rcpp::IntegerMatrix& cross_info, // columns are individuals
                             const Rcpp::NumericVector& rec_frac,   // length nrow(genotypes)-1
                             const Rcpp::IntegerVector& marker_index, // length nrow(genotypes)
                             const double error_prob,
                             const int n_draws); // number of imputations

#endif // HMM_SIMGENO_H
