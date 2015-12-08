// calculate genetic similarity (kinship matrix) from genotype probabilities
#ifndef CALC_KINSHIP_H
#define CALC_KINSHIP_H

#include <Rcpp.h>

Rcpp::NumericMatrix calc_kinship(const Rcpp::NumericVector& prob_array); // array as n_pos x n_gen x n_ind

#endif // CALC_KINSHIP_H
