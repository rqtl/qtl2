// fit a single-QTL model at a single position by Haley-Knott regression
#ifndef FIT1_HK_H
#define FIT1_HK_H

#include <Rcpp.h>

// Fit a single-QTL model at a single position
//
// genoprobs = matrix of genotype probabilities (individuals x genotypes)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list with coef, fitted, resid, rss, sigma, rank, df, SE
//
Rcpp::List fit1_hk_addcovar(const Rcpp::NumericMatrix& genoprobs,
                            const Rcpp::NumericVector& pheno,
                            const Rcpp::NumericMatrix& addcovar,
                            const Rcpp::NumericVector& weights,
                            const bool se,
                            const double tol);

// Fit a single-QTL model at a single position, with interactive covariates
//
// genoprobs = matrix of genotype probabilities (individuals x genotypes)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates
// intcovar  = interactive covariates (should also be included in addcovar)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list with coef, fitted, resid, rss, sigma, rank, df, SE
//
Rcpp::List fit1_hk_intcovar(const Rcpp::NumericVector& genoprobs,
                            const Rcpp::NumericVector& pheno,
                            const Rcpp::NumericMatrix& addcovar,
                            const Rcpp::NumericMatrix& intcovar,
                            const Rcpp::NumericVector& weights,
                            const bool se,
                            const double tol);

#endif // FIT1_HK_H
