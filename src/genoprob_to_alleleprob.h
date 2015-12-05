// convert genotype probabilities to allele probabilities
#ifndef GENOPROB_TO_ALLELEPROB_H
#define GENOPROB_TO_ALLELEPROB_H

#include <Rcpp.h>

Rcpp::NumericVector genoprob_to_alleleprob(const Rcpp::String& crosstype,
                                           const Rcpp::NumericVector& prob_array, // array as n_gen x n_ind x n_pos
                                           const bool is_x_chr);

#endif // GENOPROB_TO_ALLELEPROB_H
