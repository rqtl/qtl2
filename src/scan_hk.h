// genome scan by Haley-Knott regression
#ifndef SCAN_HK_H
#define SCAN_HK_H

#include <Rcpp.h>

// Scan a single chromosome with no additive covariates (not even intercept)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
// tol       = tolerance value for QR decomposition for linear regression
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
Rcpp::NumericMatrix scan_hk_onechr_nocovar(const Rcpp::NumericVector& genoprobs,
                                           const Rcpp::NumericMatrix& pheno,
                                           const double tol);


// Scan a single chromosome with additive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
Rcpp::NumericMatrix scan_hk_onechr(const Rcpp::NumericVector& genoprobs,
                                   const Rcpp::NumericMatrix& pheno,
                                   const Rcpp::NumericMatrix& addcovar,
                                   const double tol);

#endif // SCAN_HK_H
