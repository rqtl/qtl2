// scan chromosome by logistic regression just to get coefficients

#include "scan1coef_binary.h"
#include <Rcpp.h>

using namespace Rcpp;

#include "binreg.h"
#include "binreg_weighted.h"
#include "matrix.h"

// Scan a single chromosome to calculate coefficients, with additive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scancoef_binary_addcovar(const NumericVector& genoprobs,
                                       const NumericVector& pheno,
                                       const NumericMatrix& addcovar,
                                       const NumericVector& weights,
                                       const int maxit=100,
                                       const double tol=1e-6,
                                       const double qr_tol=1e-12,
                                       const double nu_max=30.0)
{
    const int n_ind = pheno.size();
    if(Rf_isNull(genoprobs.attr("dim")))
        throw std::invalid_argument("genoprobs should be a 3d array but has no dim attribute");
    const Dimension d = genoprobs.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("genoprobs should be a 3d array");
    const int n_pos = d[2];
    const int n_gen = d[1];
    const int n_weights = weights.size();
    const int n_addcovar = addcovar.cols();
    const int x_size = n_ind * n_gen;
    const int n_coef = n_gen + n_addcovar;

    if(n_ind != d[0])
        throw std::range_error("length(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("length(pheno) != nrow(addcovar)");
    if(n_weights > 0 && n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");

    NumericMatrix result(n_coef, n_pos);
    NumericMatrix X(n_ind, n_coef);

    // copy addcovar into matrix, at the end
    if(n_addcovar > 0)
        std::copy(addcovar.begin(), addcovar.end(), X.begin() + x_size);

    for(int pos=0, offset=0; pos<n_pos; pos++, offset += x_size) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // copy genoprobs for pos i into a matrix
        std::copy(genoprobs.begin() + offset, genoprobs.begin() + offset + x_size, X.begin());

        // do regression
        if(n_weights > 0)
            result(_,pos) = calc_coef_binreg_weighted(X, pheno, weights, maxit, tol, qr_tol, nu_max);
        else
            result(_,pos) = calc_coef_binreg(X, pheno, maxit, tol, qr_tol, nu_max);
    }

    return result;
}


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
//
// [[Rcpp::export]]
NumericMatrix scancoef_binary_intcovar(const NumericVector& genoprobs,
                                       const NumericVector& pheno,
                                       const NumericMatrix& addcovar,
                                       const NumericMatrix& intcovar,
                                       const NumericVector& weights,
                                       const int maxit=100,
                                       const double tol=1e-6,
                                       const double qr_tol=1e-12,
                                       const double nu_max=30.0)
{
    const int n_ind = pheno.size();
    if(Rf_isNull(genoprobs.attr("dim")))
        throw std::invalid_argument("genoprobs should be a 3d array but has no dim attribute");
    const Dimension d = genoprobs.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("genoprobs should be a 3d array");
    const int n_pos = d[2];
    const int n_gen = d[1];
    const int n_weights = weights.size();
    const int n_addcovar = addcovar.cols();
    const int n_intcovar = intcovar.cols();
    const int n_coef = n_gen + n_addcovar + (n_gen-1)*n_intcovar;
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");
    if(n_weights > 0 && n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");

    NumericMatrix result(n_coef, n_pos);

    for(int pos=0; pos<n_pos; pos++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // form X matrix
        NumericMatrix X = formX_intcovar(genoprobs, addcovar, intcovar, pos, false);

        // do regression
        if(n_weights > 0)
            result(_,pos) = calc_coef_binreg_weighted(X, pheno, weights, maxit, tol, qr_tol, nu_max);
        else
            result(_,pos) = calc_coef_binreg(X, pheno, maxit, tol, qr_tol, nu_max);
    }

    return result;
}


// Scan a single chromosome to calculate coefficients, with additive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list of two matrices, of coefficients and SEs (each genotypes x positions)
//
// [[Rcpp::export]]
List scancoefSE_binary_addcovar(const NumericVector& genoprobs,
                                const NumericVector& pheno,
                                const NumericMatrix& addcovar,
                                const NumericVector& weights,
                                const int maxit=100,
                                const double tol=1e-6,
                                const double qr_tol=1e-12,
                                const double nu_max=30.0)
{
    const int n_ind = pheno.size();
    if(Rf_isNull(genoprobs.attr("dim")))
        throw std::invalid_argument("genoprobs should be a 3d array but has no dim attribute");
    const Dimension d = genoprobs.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("genoprobs should be a 3d array");
    const int n_pos = d[2];
    const int n_gen = d[1];
    const int n_weights = weights.size();
    const int n_addcovar = addcovar.cols();
    const int x_size = n_ind * n_gen;
    const int n_coef = n_gen + n_addcovar;

    if(n_ind != d[0])
        throw std::range_error("length(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("length(pheno) != nrow(addcovar)");
    if(n_weights > 0 && n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");

    NumericMatrix coef(n_coef, n_pos);
    NumericMatrix se(n_coef, n_pos);
    NumericMatrix X(n_ind, n_coef);

    // copy addcovar into matrix
    if(n_addcovar > 0)
        std::copy(addcovar.begin(), addcovar.end(), X.begin() + x_size);

    for(int pos=0, offset=0; pos<n_pos; pos++, offset += x_size) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // copy genoprobs for pos i into a matrix
        std::copy(genoprobs.begin() + offset, genoprobs.begin() + offset + x_size, X.begin());

        // do regression
        List tmp;
        if(n_weights > 0)
            tmp = calc_coefSE_binreg_weighted(X, pheno, weights, maxit, tol, qr_tol, nu_max);
        else
            tmp = calc_coefSE_binreg(X, pheno, maxit, tol, qr_tol, nu_max);
        NumericVector tmpcoef = tmp[0];
        NumericVector tmpse = tmp[1];
        coef(_,pos) = tmpcoef;
        se(_,pos) = tmpse;
    }

    return List::create(Named("coef") = coef,
                        Named("SE") = se);
}


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
//
// [[Rcpp::export]]
List scancoefSE_binary_intcovar(const NumericVector& genoprobs,
                                const NumericVector& pheno,
                                const NumericMatrix& addcovar,
                                const NumericMatrix& intcovar,
                                const NumericVector& weights,
                                const int maxit=100,
                                const double tol=1e-6,
                                const double qr_tol=1e-12,
                                const double nu_max=30.0)
{
    const int n_ind = pheno.size();
    if(Rf_isNull(genoprobs.attr("dim")))
        throw std::invalid_argument("genoprobs should be a 3d array but has no dim attribute");
    const Dimension d = genoprobs.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("genoprobs should be a 3d array");
    const int n_pos = d[2];
    const int n_gen = d[1];
    const int n_weights = weights.size();
    const int n_addcovar = addcovar.cols();
    const int n_intcovar = intcovar.cols();
    const int n_coef = n_gen + n_addcovar + (n_gen-1)*n_intcovar;
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");
    if(n_weights > 0 && n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");

    NumericMatrix coef(n_coef, n_pos);
    NumericMatrix se(n_coef, n_pos);

    for(int pos=0; pos<n_pos; pos++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // form X matrix
        NumericMatrix X = formX_intcovar(genoprobs, addcovar, intcovar, pos, false);

        // do regression
        List tmp;
        if(n_weights > 0)
            tmp = calc_coefSE_binreg_weighted(X, pheno, weights, maxit, tol, qr_tol, nu_max);
        else
            tmp = calc_coefSE_binreg(X, pheno, maxit, tol, qr_tol, nu_max);
        NumericVector tmpcoef = tmp[0];
        NumericVector tmpse = tmp[1];
        coef(_,pos) = tmpcoef;
        se(_,pos) = tmpse;
    }

    return List::create(Named("coef") = coef,
                        Named("SE") = se);
}
