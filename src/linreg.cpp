// generic linear regression functions

// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

#include "linreg.h"
#include "linreg_eigen.h"

// Calculate vector of residual sum of squares (RSS) from linear regression of Y vs X
// [[Rcpp::export]]
NumericVector calc_rss_linreg(const NumericMatrix& X, const NumericMatrix& Y)
{
    return calc_mvrss_eigenqr(X, Y);
}

// Calculate matrix of residuals from linear regression of Y on X
// [[Rcpp::export]]
NumericMatrix calc_resid_linreg(const NumericMatrix& X, const NumericMatrix& Y)
{
    return calc_resid_eigenqr(X, Y);
}

// use calc_resid_linreg for a 3-dim array
// [[Rcpp::export]]
NumericVector calc_resid_linreg_3d(const NumericMatrix& X, const NumericVector& P)
{
    const unsigned int nrowx = X.rows();
    const unsigned int sizep = P.size();

    NumericMatrix pr(nrowx, sizep/nrowx);
    std::copy(P.begin(), P.end(), pr.begin()); // FIXME I shouldn't need to copy

    NumericMatrix result = calc_resid_eigenqr(X, pr);
    result.attr("dim") = P.attr("dim");

    return result;
}
