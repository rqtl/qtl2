// logistic regression
#ifndef BINREG_H
#define BINREG_H

#include <RcppEigen.h>

// logistic regression
// return just the log likelihood
double calc_ll_binreg(const Rcpp::NumericMatrix& X,
                      const Rcpp::NumericVector& y,
                      const int maxit,
                      const double tol,
                      const double qr_tol);

// logistic regression
// return just the coefficients
Rcpp::NumericVector calc_coef_binreg(const Rcpp::NumericMatrix& X,
                                     const Rcpp::NumericVector& y,
                                     const int maxit,
                                     const double tol,
                                     const double qr_tol);
#endif // BINREG_H
