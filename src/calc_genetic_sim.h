// calculate genetic similarity from genotype probabilities
#ifndef CALC_GENETIC_SIM_H
#define CALC_GENETIC_SIM_H

#include <Rcpp.h>

Rcpp::NumericMatrix calc_genetic_sim(const Rcpp::NumericVector& prob_array); // array as n_pos x n_gen x n_ind

#endif // CALC_GENETIC_SIM_H
