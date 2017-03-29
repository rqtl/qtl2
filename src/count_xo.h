// count number of crossovers
#ifndef COUNT_XO_H
#define COUNT_XO_H

#include <Rcpp.h>

Rcpp::IntegerVector count_xo(const Rcpp::IntegerMatrix geno, // genotype matrix markers x individuals
                             const Rcpp::String& crosstype,
                             const bool is_X_chr);

#endif // COUNT_XO_H
