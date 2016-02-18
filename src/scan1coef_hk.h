// scan chromosome by Haley-Knott regression just to get coefficients
#ifndef SCAN1COEF_HK_H
#define SCAN1COEF_HK_H

#include <Rcpp.h>

// Scan a single chromosome to calculate coefficients, with no covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
Rcpp::NumericMatrix scancoef_hk_nocovar(const Rcpp::NumericVector& genoprobs,
                                        const Rcpp::NumericVector& pheno,
                                        const Rcpp::NumericVector& weights,
                                        const double tol);

// Scan a single chromosome to calculate coefficients, with additive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
Rcpp::NumericMatrix scancoef_hk_addcovar(const Rcpp::NumericVector& genoprobs,
                                         const Rcpp::NumericVector& pheno,
                                         const Rcpp::NumericMatrix& addcovar,
                                         const Rcpp::NumericVector& weights,
                                         const double tol);

// Scan a single chromosome to calculate coefficients, with interactive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// intcovar  = interactive covariates (should also be included in addcovar)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
Rcpp::NumericMatrix scancoef_hk_intcovar(const Rcpp::NumericVector& genoprobs,
                                         const Rcpp::NumericVector& pheno,
                                         const Rcpp::NumericMatrix& addcovar,
                                         const Rcpp::NumericMatrix& intcovar,
                                         const Rcpp::NumericVector& weights,
                                         const double tol);




// Scan a single chromosome to calculate coefficients, with no covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
Rcpp::List scancoefSE_hk_nocovar(const Rcpp::NumericVector& genoprobs,
                                 const Rcpp::NumericVector& pheno,
                                 const Rcpp::NumericVector& weights,
                                 const double tol);


// Scan a single chromosome to calculate coefficients, with additive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
Rcpp::List scancoefSE_hk_addcovar(const Rcpp::NumericVector& genoprobs,
                                  const Rcpp::NumericVector& pheno,
                                  const Rcpp::NumericMatrix& addcovar,
                                  const Rcpp::NumericVector& weights,
                                  const double tol);


// Scan a single chromosome to calculate coefficients, with interactive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// intcovar  = interactive covariates (should also be included in addcovar)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
Rcpp::List scancoefSE_hk_intcovar(const Rcpp::NumericVector& genoprobs,
                                  const Rcpp::NumericVector& pheno,
                                  const Rcpp::NumericMatrix& addcovar,
                                  const Rcpp::NumericMatrix& intcovar,
                                  const Rcpp::NumericVector& weights,
                                  const double tol);
#endif // SCAN1COEF_HK_H
