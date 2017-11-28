// locate crossovers
#ifndef LOCATE_XO_H
#define LOCATE_XO_H

#include <Rcpp.h>

Rcpp::List locate_xo(const Rcpp::IntegerMatrix geno, // genotype matrix markers x individuals
                     const Rcpp::NumericVector map,
                     const Rcpp::String& crosstype,
                     const bool is_X_chr);

#endif // LOCATE_XO_H
