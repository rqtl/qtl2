// interpolate genotype probabilities
#ifndef INTERP_GENOPROB_H
#define INTERP_GENOPROB_H

#include <Rcpp.h>

Rcpp::NumericVector interp_genoprob_onechr(const Rcpp::NumericVector& genoprob,
                                           const Rcpp::NumericVector& map,
                                           const Rcpp::IntegerVector& pos_index);

#endif // INTERP_GENOPROB_H
