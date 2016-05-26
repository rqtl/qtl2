// identify genotype with maximum marginal probability, arg max Pr(g)

#include <Rcpp.h>

Rcpp::IntegerMatrix maxmarg(const Rcpp::NumericVector& prob_array, const int minprob);
