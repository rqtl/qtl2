// scan chromosome by logistic regression just to get coefficients
#ifndef SCAN1COEF_BINARY_H
#define SCAN1COEF_BINARY_H

#include <Rcpp.h>

// Scan a single chromosome to calculate coefficients, with additive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
Rcpp::NumericMatrix scancoef_binary_addcovar(const Rcpp::NumericVector& genoprobs,
                                             const Rcpp::NumericVector& pheno,
                                             const Rcpp::NumericMatrix& addcovar,
                                             const Rcpp::NumericVector& weights,
                                             const int maxit,
                                             const double tol,
                                             const double qr_tol,
                                             const double nu_max);

// Scan a single chromosome to calculate coefficients, with interactive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates
// intcovar  = interactive covariates (should also be included in addcovar)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
Rcpp::NumericMatrix scancoef_binary_intcovar(const Rcpp::NumericVector& genoprobs,
                                             const Rcpp::NumericVector& pheno,
                                             const Rcpp::NumericMatrix& addcovar,
                                             const Rcpp::NumericMatrix& intcovar,
                                             const Rcpp::NumericVector& weights,
                                             const int maxit,
                                             const double tol,
                                             const double qr_tol,
                                             const double nu_max);

// Scan a single chromosome to calculate coefficients, with additive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list of two matrices, of coefficients and SEs (each genotypes x positions)
Rcpp::List scancoefSE_binary_addcovar(const Rcpp::NumericVector& genoprobs,
                                      const Rcpp::NumericVector& pheno,
                                      const Rcpp::NumericMatrix& addcovar,
                                      const Rcpp::NumericVector& weights,
                                      const int maxit,
                                      const double tol,
                                      const double qr_tol,
                                      const double nu_max);


// Scan a single chromosome to calculate coefficients, with interactive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates
// intcovar  = interactive covariates (should also be included in addcovar)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list of two matrices, of coefficients and SEs (each genotypes x positions)
Rcpp::List scancoefSE_binary_intcovar(const Rcpp::NumericVector& genoprobs,
                                      const Rcpp::NumericVector& pheno,
                                      const Rcpp::NumericMatrix& addcovar,
                                      const Rcpp::NumericMatrix& intcovar,
                                      const Rcpp::NumericVector& weights,
                                      const int maxit,
                                      const double tol,
                                      const double qr_tol,
                                      const double nu_max);

#endif // SCAN1COEF_BINARY_H
