// fit a single-QTL model at a single position by LMM
#ifndef FIT1_PG_H
#define FIT1_PG_H

#include <Rcpp.h>


// fit single-QTL model at a single position
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// eigenvec  = eigenvectors from eigen decomposition of kinship matrix
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list with a bunch of stuff
//
Rcpp::List fit1_pg_addcovar(const Rcpp::NumericMatrix& genoprobs,
                            const Rcpp::NumericVector& pheno,
                            const Rcpp::NumericMatrix& addcovar,
                            const Rcpp::NumericMatrix& eigenvec,
                            const Rcpp::NumericVector& weights,
                            const bool se,
                            const double tol);


// fit single-QTL model at a single position
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// intcovar  = interactive covariates (should also be included in addcovar)
// eigenvec  = eigenvectors from eigen decomposition of kinship matrix
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list of a bunch of stuff
//
Rcpp::List fit1_pg_intcovar(const Rcpp::NumericMatrix& genoprobs,
                            const Rcpp::NumericVector& pheno,
                            const Rcpp::NumericMatrix& addcovar,
                            const Rcpp::NumericMatrix& intcovar,
                            const Rcpp::NumericMatrix& eigenvec,
                            const Rcpp::NumericVector& weights,
                            const bool se,
                            const double tol);

#endif // FIT1_PG_H
