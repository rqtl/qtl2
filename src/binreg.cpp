// logistic regression via RcppEigen

// [[Rcpp::depends(RcppEigen)]]

#include "binreg.h"
#include <RcppEigen.h>
#include "binreg_eigen.h"

using namespace Rcpp;
using namespace Eigen;

// logistic regression
// return just the log likelihood
// [[Rcpp::export]]
double calc_ll_binreg(const NumericMatrix& X, const NumericVector& y,
                      const int maxit=100, const double tol=1e-6,
                      const double qr_tol=1e-12, const double eta_max=30.0)
{
    return calc_ll_binreg_eigenqr(X, y, maxit, tol, qr_tol, eta_max);
}

// logistic regression
// return just the coefficients
// [[Rcpp::export]]
NumericVector calc_coef_binreg(const NumericMatrix& X, const NumericVector& y,
                               const int maxit=100, const double tol=1e-6,
                               const double qr_tol=1e-12, const double eta_max=30.0)
{
    return calc_coef_binreg_eigenqr(X, y, maxit, tol, qr_tol, eta_max);
}

// logistic regression
// return the coefficients and SEs
// [[Rcpp::export]]
List calc_coefSE_binreg(const NumericMatrix& X, const NumericVector& y,
                        const int maxit=100, const double tol=1e-6,
                        const double qr_tol=1e-12, const double eta_max=30.0)
{
    return calc_coefSE_binreg_eigenqr(X, y, maxit, tol, qr_tol, eta_max);
}

// logistic regression
// return (llik, individual contributions to llik, fitted probabilities, coef, SE
// [[Rcpp::export]]
List fit_binreg(const NumericMatrix& X, const NumericVector& y,
                const bool se=true, // whether to include SEs
                const int maxit=100, const double tol=1e-6,
                const double qr_tol=1e-12, const double eta_max=30.0)
{
    return fit_binreg_eigenqr(X, y, se, maxit, tol, qr_tol, eta_max);
}
