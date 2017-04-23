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
                      const double qr_tol=1e-12)
{
    return calc_ll_binreg_eigenqr(X, y, maxit, tol, qr_tol);
}
