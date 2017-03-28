// calculate genotyping error LOD scores
// (assumes constant is_female and cross_info and pre-calcs the step and emit matrices)
#ifndef HMM_CALCERRORLOD_H
#define HMM_CALCERRORLOD_H

#include <Rcpp.h>

Rcpp::NumericMatrix calc_errorlod(const Rcpp::String& crosstype,
                                  const Rcpp::NumericVector& probs, // genotype probs [genotype, ind, marker]
                                  const Rcpp::IntegerMatrix& genotypes, // columns are individuals, rows are markers
                                  const Rcpp::IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                                  const bool is_X_chr,
                                  const bool is_female, // same for all individuals
                                  const Rcpp::IntegerVector& cross_info); // same for all individuals

#endif // HMM_CALCERRORLOD_H
