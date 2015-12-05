// forward-backward equations for HMM
#ifndef HMM_FORWBACK_H
#define HMM_FORWBACK_H

#include <Rcpp.h>
#include "cross.h"

// forward equations
Rcpp::NumericMatrix forwardEquations(QTLCross* cross,
                                     const Rcpp::IntegerVector& genotypes,
                                     const Rcpp::IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                                     const bool is_X_chr,
                                     const bool is_female,
                                     const Rcpp::IntegerVector& cross_info,
                                     const Rcpp::NumericVector& rec_frac,
                                     const Rcpp::IntegerVector& marker_index,
                                     const double error_prob,
                                     const Rcpp::IntegerVector& poss_gen);


// backward Equations
Rcpp::NumericMatrix backwardEquations(QTLCross* cross,
                                      const Rcpp::IntegerVector& genotypes,
                                      const Rcpp::IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                                      const bool is_X_chr,
                                      const bool is_female,
                                      const Rcpp::IntegerVector& cross_info,
                                      const Rcpp::NumericVector& rec_frac,
                                      const Rcpp::IntegerVector& marker_index,
                                      const double error_prob,
                                      const Rcpp::IntegerVector& poss_gen);

#endif // HMM_FORWBACK_H
