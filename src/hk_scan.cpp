// simple genome scan by linear regression

#include <Rcpp.h>

using namespace Rcpp;

#include "linreg.h"
#include "hk_scan.h"

// use calc_resid_lapack for a 3-dim array
// [[Rcpp::export]]
NumericVector calc_resid_linreg_3d(const NumericMatrix& X, const NumericVector& P)
{
    int nrowx = X.rows();
    int sizep = P.size();

    NumericMatrix pr(nrowx, sizep/nrowx);
    std::copy(P.begin(), P.end(), pr.begin());

    NumericMatrix result = calc_resid_linreg(X, pr);
    result.attr("dim") = P.attr("dim");

    return result;
}
// .calc_resid(X, aperm(pr[[1]], c(1,3,2))) // then permute genotypes again
