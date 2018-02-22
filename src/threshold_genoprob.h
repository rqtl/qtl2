// threshold genoprobs, setting small values to 0
#ifndef THRESHOLD_GENOPROB_H
#define THRESHOLD_GENOPROB_H

#include <Rcpp.h>

// threshold genoprobs, setting small values to 0
Rcpp::NumericVector threshold_genoprob(const Rcpp::NumericVector& prob_array, // array as n_gen x n_ind x n_pos
                                       const double threshold);

#endif // THRESHOLD_GENOPROB_H
