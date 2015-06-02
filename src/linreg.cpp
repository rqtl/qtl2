// generic linear regression functions

#include <Rcpp.h>

using namespace Rcpp;

#include "linreg.h"
#include "linreg_lapack.h"

// Calculate vector of residual sum of squares (RSS) from linear regression of Y vs X
// [[Rcpp::export]]
NumericVector calc_rss_linreg(const NumericMatrix& X, const NumericMatrix& Y)
{
    return calc_rss_lapack(X, Y,
                           false, // skip_dgels
                           1e-10); // tolerance
}

// Calculate matrix of residuals from linear regression of Y on X
// [[Rcpp::export]]
NumericMatrix calc_resid_linreg(const NumericMatrix& X, const NumericMatrix& Y)
{
    return calc_resid_lapack(X, Y,
                             false, // skip dgels
                             1e-10); // tolerance
}

// use calc_resid_linreg for a 3-dim array
// [[Rcpp::export]]
NumericVector calc_resid_linreg_3d(const NumericMatrix& X, const NumericVector& P)
{
    int nrowx = X.rows();
    int sizep = P.size();

    NumericMatrix pr(nrowx, sizep/nrowx);
    std::copy(P.begin(), P.end(), pr.begin()); // FIXME I shouldn't need to copy

    NumericMatrix result = calc_resid_linreg(X, pr);
    result.attr("dim") = P.attr("dim");

    return result;
}
// calc_resid_linreg_3d(X, aperm(pr[[1]], c(1,3,2))) // then permute genotypes again
