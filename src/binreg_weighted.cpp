// (weighted) logistic regression via RcppEigen

// [[Rcpp::depends(RcppEigen)]]

#include "binreg_weighted.h"
#include <RcppEigen.h>
#include "binreg_weighted_eigen.h"

using namespace Rcpp;
using namespace Eigen;

// logistic regression
// return just the log likelihood
// this version uses weights
// [[Rcpp::export]]
double calc_ll_binreg_weighted(const NumericMatrix& X, const NumericVector& y,
                               const NumericVector &weights,
                               const int maxit=100, const double tol=1e-6,
                               const double qr_tol=1e-12, const double nu_max=30.0)
{
    return calc_ll_binreg_weighted_eigenqr(X, y, weights, maxit, tol, qr_tol, nu_max);
}

// logistic regression
// return just the coefficients
// this version uses weights
// [[Rcpp::export]]
NumericVector calc_coef_binreg_weighted(const NumericMatrix& X, const NumericVector& y,
                                        const NumericVector &weights,
                                        const int maxit=100, const double tol=1e-6,
                                        const double qr_tol=1e-12, const double nu_max=30.0)
{
    return calc_coef_binreg_weighted_eigenqr(X, y, weights, maxit, tol, qr_tol, nu_max);
}

// logistic regression
// return the coefficients and SEs
// this version uses weights
// [[Rcpp::export]]
List calc_coefSE_binreg_weighted(const NumericMatrix& X, const NumericVector& y,
                                 const NumericVector &weights,
                                 const int maxit=100, const double tol=1e-6,
                                 const double qr_tol=1e-12, const double nu_max=30.0)
{
    return calc_coefSE_binreg_weighted_eigenqr(X, y, weights, maxit, tol, qr_tol, nu_max);
}

// logistic regression
// return (llik, individual contributions to llik, fitted probabilities, coef, SE
// this version uses weights
// [[Rcpp::export]]
List fit_binreg_weighted(const NumericMatrix& X, const NumericVector& y,
                         const NumericVector &weights,
                         const bool se=true, // whether to include SEs
                         const int maxit=100, const double tol=1e-6,
                         const double qr_tol=1e-12, const double nu_max=30.0)
{
    return fit_binreg_weighted_eigenqr(X, y, weights, se, maxit, tol, qr_tol, nu_max);
}
