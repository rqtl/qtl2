// scan chromosome by Haley-Knott regression just to get coefficients

#include "scan1coef_hk.h"
#include <Rcpp.h>

using namespace Rcpp;

#include "linreg.h"
#include "matrix.h"

// Scan a single chromosome to calculate coefficients, with no covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scancoef_hk_nocovar(const NumericVector& genoprobs,
                                  const NumericVector& pheno,
                                  const NumericVector& weights,
                                  const double tol=1e-12)
{
    const unsigned int n_ind = pheno.size();
    const Dimension d = genoprobs.attr("dim");
    const unsigned int n_pos = d[2];
    const unsigned int n_gen = d[1];
    const unsigned int x_size = n_ind * n_gen;
    const unsigned int n_weights = weights.size();

    if(n_ind != d[0])
        throw std::range_error("length(pheno) != nrow(genoprobs)");
    if(n_weights > 0 && n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");

    NumericMatrix result(n_gen, n_pos);
    NumericMatrix X(n_ind, n_gen);

    for(unsigned int pos=0, offset=0; pos<n_pos; pos++, offset += x_size) {
        // copy genoprobs for pos i into a matrix
        std::copy(genoprobs.begin() + offset, genoprobs.begin() + offset + x_size, X.begin());

        // multiply by square-root weights, if necessary
        if(n_weights > 0) X = weighted_matrix(X, weights);

        // do regression
        result(_,pos) = calc_coef_linreg(X, pheno, tol);
    }

    return result;
}


// Scan a single chromosome to calculate coefficients, with additive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scancoef_hk_addcovar(const NumericVector& genoprobs,
                                   const NumericVector& pheno,
                                   const NumericMatrix& addcovar,
                                   const NumericVector& weights,
                                   const double tol=1e-12)
{
    const unsigned int n_ind = pheno.size();
    const Dimension d = genoprobs.attr("dim");
    const unsigned int n_pos = d[2];
    const unsigned int n_gen = d[1];
    const unsigned int n_weights = weights.size();
    const unsigned int n_addcovar = addcovar.cols();
    const unsigned int x_size = n_ind * n_gen;
    const unsigned int n_coef = n_gen + n_addcovar;

    if(n_ind != d[0])
        throw std::range_error("length(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("length(pheno) != nrow(addcovar)");
    if(n_weights > 0 && n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");

    NumericMatrix result(n_coef, n_pos);
    NumericMatrix X(n_ind, n_coef);

    for(unsigned int pos=0, offset=0; pos<n_pos; pos++, offset += x_size) {
        // copy genoprobs for pos i into a matrix
        std::copy(genoprobs.begin() + offset, genoprobs.begin() + offset + x_size, X.begin());

        // copy addcovar into matrix
        std::copy(addcovar.begin(), addcovar.end(), X.begin() + x_size);

        // multiply by square-root weights, if necessary
        if(n_weights > 0) X = weighted_matrix(X, weights);

        // do regression
        result(_,pos) = calc_coef_linreg(X, pheno, tol);
    }

    return result;
}


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
//
// [[Rcpp::export]]
NumericMatrix scancoef_hk_intcovar(const NumericVector& genoprobs,
                                   const NumericVector& pheno,
                                   const NumericMatrix& addcovar,
                                   const NumericMatrix& intcovar,
                                   const NumericVector& weights,
                                   const double tol=1e-12)
{
    const unsigned int n_ind = pheno.size();
    const Dimension d = genoprobs.attr("dim");
    const unsigned int n_pos = d[2];
    const unsigned int n_gen = d[1];
    const unsigned int n_weights = weights.size();
    const unsigned int n_addcovar = addcovar.cols();
    const unsigned int n_intcovar = intcovar.cols();
    const unsigned int n_coef = n_gen + n_addcovar + (n_gen-1)*n_intcovar;
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");
    if(n_weights > 0 && n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");

    NumericMatrix result(n_coef, n_pos);

    for(unsigned int pos=0; pos<n_pos; pos++) {
        // form X matrix
        NumericMatrix X = formX_intcovar(genoprobs, addcovar, intcovar, pos);
        if(n_weights > 0) X = weighted_matrix(X, weights);

        // do regression
        result(_,pos) = calc_coef_linreg(X, pheno, tol);
    }

    return result;
}
