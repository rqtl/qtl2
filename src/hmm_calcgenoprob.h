// calculate conditional genotype probabilities given multipoint marker data
#ifndef HMM_CALCGENOPROB_H
#define HMM_CALCGENOPROB_H

#include <Rcpp.h>

Rcpp::NumericVector calc_genoprob(const Rcpp::String& crosstype,
                                  const Rcpp::IntegerMatrix& genotypes, // columns are individuals, rows are markers
                                  const Rcpp::IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                                  const bool is_X_chr,
                                  const Rcpp::LogicalVector& is_female, // length n_ind
                                  const Rcpp::IntegerMatrix& cross_info, // columns are individuals
                                  const Rcpp::NumericVector& rec_frac,   // length nrow(genotypes)-1
                                  const Rcpp::IntegerVector& marker_index, // length nrow(genotypes)
                                  const double error_prob);

#endif // HMM_CALCGENOPROB_H
