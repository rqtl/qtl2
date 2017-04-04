// fit a single-QTL model at a single position by Haley-Knott regression

#include "fit1_hk.h"
#include <Rcpp.h>

using namespace Rcpp;

#include "linreg.h"
#include "matrix.h"

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
// [[Rcpp::export]]
List fit1_hk_addcovar(const NumericMatrix& genoprobs,
                      const NumericVector& pheno,
                      const NumericMatrix& addcovar,
                      const NumericVector& weights,
                      const bool se,
                      const double tol=1e-12)
{
    const int n_ind = pheno.size();
    const int n_gen = genoprobs.cols();
    const int n_weights = weights.size();
    const int n_addcovar = addcovar.cols();
    const int x_size = n_ind * n_gen;
    const int n_coef = n_gen + n_addcovar;

    if(n_ind != genoprobs.rows())
        throw std::range_error("length(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("length(pheno) != nrow(addcovar)");
    if(n_weights > 0 && n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");

    NumericMatrix X(n_ind, n_coef);

    // copy genoprobs into matrix
    std::copy(genoprobs.begin(), genoprobs.end(), X.begin());

    // copy addcovar into matrix
    if(n_addcovar > 0)
        std::copy(addcovar.begin(), addcovar.end(), X.begin() + x_size);

    // multiply by square-root weights, if necessary
    if(n_weights > 0) X = weighted_matrix(X, weights);

    return fit_linreg(X, pheno, se, tol);
}


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
// [[Rcpp::export]]
List fit1_hk_intcovar(const NumericMatrix& genoprobs,
                      const NumericVector& pheno,
                      const NumericMatrix& addcovar,
                      const NumericMatrix& intcovar,
                      const NumericVector& weights,
                      const bool se,
                      const double tol=1e-12)
{
    const int n_ind = pheno.size();
    const int n_weights = weights.size();

    if(n_ind != genoprobs.rows())
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");
    if(n_weights > 0 && n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");

    // form X matrix
    NumericMatrix X = formX_intcovar(genoprobs, addcovar, intcovar, 0, false);
    if(n_weights > 0) X = weighted_matrix(X, weights);

    return fit_linreg(X, pheno, se, tol);
}
