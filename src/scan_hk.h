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

// Scan a single chromosome with additive covariates and weights
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of (weighted) residual sums of squares (RSS) (phenotypes x positions)
Rcpp::NumericMatrix scan_hk_onechr_weighted(const Rcpp::NumericVector& genoprobs,
                                            const Rcpp::NumericMatrix& pheno,
                                            const Rcpp::NumericMatrix& addcovar,
                                            const Rcpp::NumericVector& weights,
                                            const double tol);

// Scan a single chromosome with interactive covariates
// this version should be fast but requires more memory
// (since we first expand the genotype probabilities to probs x intcovar)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
// intcovar  = interactive covariates (should also be included in addcovar)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
Rcpp::NumericMatrix scan_hk_onechr_intcovar_highmem(const Rcpp::NumericVector& genoprobs,
                                                    const Rcpp::NumericMatrix& pheno,
                                                    const Rcpp::NumericMatrix& addcovar,
                                                    const Rcpp::NumericMatrix& intcovar,
                                                    const double);

#endif // SCAN_HK_H
