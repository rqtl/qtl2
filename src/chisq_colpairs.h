// perform chi-square tests on all pairs of columns of a matrix
#ifndef CHISQ_COLPAIRS_H
#define CHISQ_COLPAIRS_H

#include <Rcpp.h>

Rcpp::NumericMatrix chisq_colpairs(const Rcpp::IntegerMatrix& matrix); // matrix of integers; should be contiguous

#endif // CHISQ_COLPAIRS_H
