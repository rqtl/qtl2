// logistic regression via RcppEigen
#ifndef BINREG_EIGEN_H
#define BINREG_EIGEN_H

#include <RcppEigen.h>

// logistic regression by "LLt" Cholesky decomposition
// return just the log likelihood
double calc_ll_binreg_eigenchol(const Rcpp::NumericMatrix& X,
                                const Rcpp::NumericVector& y,
                                const int maxit,
                                const double tol);

// logistic regression by Qr decomposition with column pivoting
// return just the log likelihood
double calc_ll_binreg_eigenqr(const Rcpp::NumericMatrix& X,
                              const Rcpp::NumericVector& y,
                              const int maxit,
                              const double tol,
                              const double qr_tol);

// logistic regression
// return just the coefficients
Rcpp::NumericVector calc_coef_binreg_eigenqr(const Rcpp::NumericMatrix& X,
                                             const Rcpp::NumericVector& y,
                                             const int maxit,
                                             const double tol,
                                             const double qr_tol);

#endif // BINREG_EIGEN_H
