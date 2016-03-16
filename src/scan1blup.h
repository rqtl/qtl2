// scan chromosome to get BLUPs of coefficients
#ifndef SCAN1BLUP_H
#define SCAN1BLUP_H

#include <Rcpp.h>

// Scan a single chromosome to get BLUPs of coefficients
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (must include intercept)
// se        = If TRUE, calculate SEs
// reml      = If TRUE, use REML to estimate variance components; otherwise use maximum
//             likelihood
// preserve_intercept = If FALSE, add the intercept to the BLUPs and remove that column
// tol       = Numeric tolerance
//
// output    = List with two matrices, of coefficients and SEs (each coefficients x positions)
Rcpp::List scanblup(const Rcpp::NumericVector& genoprobs,
                    const Rcpp::NumericVector& pheno,
                    const Rcpp::NumericMatrix& addcovar,
                    const bool se,
                    const bool reml,
                    const bool preserve_intercept,
                    const double tol);

#endif // SCAN1BLUP_H
