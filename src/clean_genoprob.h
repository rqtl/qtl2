// clean genoprobs, setting small values to 0
#ifndef CLEAN_GENOPROB_H
#define CLEAN_GENOPROB_H

#include <Rcpp.h>

// clean genoprobs, setting small values to 0
Rcpp::NumericVector clean_genoprob(const Rcpp::NumericVector& prob_array, // array as n_ind x n_gen x n_pos
                                   double value_threshold,
                                   double column_threshold);

#endif // CLEAN_GENOPROB_H
