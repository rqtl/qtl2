// fit a single-QTL model at a single position by Haley-Knott regression
#ifndef FIT1_BINARY_H
#define FIT1_BINARY_H

#include <Rcpp.h>

// Fit a single-QTL model at a single position
//
// genoprobs = matrix of genotype probabilities (individuals x genotypes)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed; values in [0,1])
// addcovar  = additive covariates
// weights   = vector of weights
//
// output    = list with lod, fitted probabilities, coef, SE
Rcpp::List fit1_binary_addcovar(const Rcpp::NumericMatrix& genoprobs,
                                const Rcpp::NumericVector& pheno,
                                const Rcpp::NumericMatrix& addcovar,
                                const Rcpp::NumericVector& weights,
                                const bool se,
                                const int maxit,
                                const double tol,
                                const double qr_tol);


// Fit a single-QTL model at a single position, with interactive covariates
//
// genoprobs = matrix of genotype probabilities (individuals x genotypes)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed; values in [0,1])
// addcovar  = additive covariates
// intcovar  = interactive covariates (should also be included in addcovar)
// weights   = vector of weights
Rcpp::List fit1_binary_intcovar(const Rcpp::NumericMatrix& genoprobs,
                                const Rcpp::NumericVector& pheno,
                                const Rcpp::NumericMatrix& addcovar,
                                const Rcpp::NumericMatrix& intcovar,
                                const Rcpp::NumericVector& weights,
                                const bool se,
                                const int maxit,
                                const double tol,
                                const double qr_tol);

#endif // FIT1_BINARY_H
