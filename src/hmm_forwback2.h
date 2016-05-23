// forward-backward equations for HMM
// (this version assumes constant is_female and cross_info and pre-calcs the step and emit matrices)
#ifndef HMM_FORWBACK2_H
#define HMM_FORWBACK2_H

#include <Rcpp.h>
#include "cross.h"

// forward equations
Rcpp::NumericMatrix forwardEquations2(const Rcpp::IntegerVector& genotypes,
                                      const Rcpp::NumericVector& init_vector,
                                      const std::vector<Rcpp::NumericMatrix>& emit_matrix,
                                      const std::vector<Rcpp::NumericMatrix>& step_matrix,
                                      const Rcpp::IntegerVector& marker_index,
                                      const Rcpp::IntegerVector& poss_gen);


// backward Equations
Rcpp::NumericMatrix backwardEquations2(const Rcpp::IntegerVector& genotypes,
                                       const Rcpp::NumericVector& init_vector,
                                      const std::vector<Rcpp::NumericMatrix>& emit_matrix,
                                      const std::vector<Rcpp::NumericMatrix>& step_matrix,
                                      const Rcpp::IntegerVector& marker_index,
                                      const Rcpp::IntegerVector& poss_gen);

#endif // HMM_FORWBACK2_H
