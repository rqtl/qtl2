// fit a single-QTL model at a single position by Haley-Knott regression

#include "fit1_binary.h"
#include <Rcpp.h>

using namespace Rcpp;

#include "binreg.h"
#include "binreg_weighted.h"
#include "matrix.h"

// Fit a single-QTL model at a single position
//
// genoprobs = matrix of genotype probabilities (individuals x genotypes)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed; values in [0,1])
// addcovar  = additive covariates
// weights   = vector of weights
//
// output    = list with lod, fitted probabilities, coef, SE
//
// [[Rcpp::export]]
List fit1_binary_addcovar(const NumericMatrix& genoprobs,
                          const NumericVector& pheno,
                          const NumericMatrix& addcovar,
                          const NumericVector& weights,
                          const bool se=false,
                          const int maxit=100,
                          const double tol=1e-6,
                          const double qr_tol=1e-12,
                          const double nu_max=30.0)
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

    if(n_weights > 0)
        return fit_binreg_weighted(X, pheno, weights, se, maxit, tol, qr_tol, nu_max);
    else
        return fit_binreg(X, pheno, se, maxit, tol, qr_tol, nu_max);
}


// Fit a single-QTL model at a single position, with interactive covariates
//
// genoprobs = matrix of genotype probabilities (individuals x genotypes)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed; values in [0,1])
// addcovar  = additive covariates
// intcovar  = interactive covariates (should also be included in addcovar)
// weights   = vector of weights
//
// output    = list with coef, fitted, resid, rss, sigma, rank, df, SE
//
// [[Rcpp::export]]
List fit1_binary_intcovar(const NumericMatrix& genoprobs,
                          const NumericVector& pheno,
                          const NumericMatrix& addcovar,
                          const NumericMatrix& intcovar,
                          const NumericVector& weights,
                          const bool se=true,
                          const int maxit=100,
                          const double tol=1e-6,
                          const double qr_tol=1e-12,
                          const double nu_max=30.0)
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

    if(n_weights > 0)
        return fit_binreg_weighted(X, pheno, weights, se, maxit, tol, qr_tol, nu_max);
    else
        return fit_binreg(X, pheno, se, maxit, tol, qr_tol, nu_max);
}
