// generic linear regression functions

// [[Rcpp::depends(RcppEigen)]]

#include "linreg.h"
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

#include "linreg_eigen.h"

// Calculate vector of residual sum of squares (RSS) from linear regression of Y on X
// [[Rcpp::export]]
NumericVector calc_rss_linreg(const NumericMatrix& X, const NumericMatrix& Y,
                              const double tol=1e-12)
{
    return calc_mvrss_eigenqr(X, Y, tol);
}

// Calculate just the coefficients from linear regression of y on X
// [[Rcpp::export]]
NumericVector calc_coef_linreg(const NumericMatrix& X, const NumericVector& y,
                               const double tol=1e-12)
{
    return calc_coef_linreg_eigenqr(X, y, tol);
}

// Calculate coefficients and SEs from linear regression of y on X
// [[Rcpp::export]]
List calc_coefSE_linreg(const NumericMatrix& X, const NumericVector& y,
                        const double tol=1e-12)
{
    return calc_coefSE_linreg_eigenqr(X, y, tol);
}

// Calculate matrix of residuals from linear regression of Y on X
// [[Rcpp::export]]
NumericMatrix calc_resid_linreg(const NumericMatrix& X, const NumericMatrix& Y,
                                const double tol=1e-12)
{
    return calc_resid_eigenqr(X, Y, tol);
}

// use calc_resid_linreg for a 3-dim array
// [[Rcpp::export]]
NumericVector calc_resid_linreg_3d(const NumericMatrix& X, const NumericVector& P,
                                   const double tol=1e-12)
{
    const int nrowx = X.rows();
    if(Rf_isNull(P.attr("dim")))
        throw std::invalid_argument("P should be a 3d array but has no dim attribute");
    const Dimension d = P.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("P should be a 3d array");
    if(d[0] != nrowx)
        throw std::range_error("nrow(X) != nrow(P)");

    NumericMatrix pr(nrowx, d[1]*d[2]);
    std::copy(P.begin(), P.end(), pr.begin()); // FIXME I shouldn't need to copy

    NumericMatrix result = calc_resid_eigenqr(X, pr, tol);
    result.attr("dim") = d;

    return result;
}

// least squares, returning everything
// output is list of (coef, fitted, resid, rss, sigma, rank, df, SE)
//
// argument se indicates whether to calculate standard errors (SE)
//
// [[Rcpp::export]]
List fit_linreg(const NumericMatrix& X, const NumericVector& y,
                const bool se=true, const double tol=1e-12)
{
    return fit_linreg_eigenqr(X, y, se, tol);
}
