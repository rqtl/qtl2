// scan chromosome with linear mixed model
#ifndef SCAN_LMM_H
#define SCAN_LMM_H

#include <Rcpp.h>

// REML scan of a single chromosome with additive covariates and weights
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix with one column of numeric phenotypes
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of log restricted likelihood values, but off by -sum(log(weights))/2
Rcpp::NumericVector scan_reml_onechr(const Rcpp::NumericVector& genoprobs,
                                     const Rcpp::NumericMatrix& pheno,
                                     const Rcpp::NumericMatrix& addcovar,
                                     const Rcpp::NumericVector& weights,
                                     const double tol);

// REML scan of a single chromosome with interactive covariates
// this version should be fast but requires more memory
// (since we first expand the genotype probabilities to probs x intcovar)
// and this one allows weights for the individuals (the same for all phenotypes)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix with one column of numeric phenotypes
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
// intcovar  = interactive covariates (should also be included in addcovar)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
Rcpp::NumericVector scan_reml_onechr_intcovar_highmem(const Rcpp::NumericVector& genoprobs,
                                                      const Rcpp::NumericMatrix& pheno,
                                                      const Rcpp::NumericMatrix& addcovar,
                                                      const Rcpp::NumericMatrix& intcovar,
                                                      const Rcpp::NumericVector& weights,
                                                      const double tol);

// REML scan of a single chromosome with interactive covariates
// this version uses less memory but will be slower
// (since we need to work with each position, one at a time)
// and this one allows weights for the individuals (the same for all phenotypes)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix with one column of numeric phenotypes
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
// intcovar  = interactive covariates (should also be included in addcovar)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
Rcpp::NumericVector scan_reml_onechr_intcovar_lowmem(const Rcpp::NumericVector& genoprobs,
                                                     const Rcpp::NumericMatrix& pheno,
                                                     const Rcpp::NumericMatrix& addcovar,
                                                     const Rcpp::NumericMatrix& intcovar,
                                                     const Rcpp::NumericVector& weights,
                                                     const double tol);

// calculate logdetXpX many times, along positions of genotype array
Rcpp::NumericVector calc_logdetXpX_many(const Rcpp::NumericVector& genoprobs,
                                        const Rcpp::NumericMatrix& addcovar);

#endif // SCAN_LMM_H
