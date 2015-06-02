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
