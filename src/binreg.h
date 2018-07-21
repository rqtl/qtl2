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
                      const double qr_tol,
                      const double nu_max);

// logistic regression
// return just the coefficients
Rcpp::NumericVector calc_coef_binreg(const Rcpp::NumericMatrix& X,
                                     const Rcpp::NumericVector& y,
                                     const int maxit,
                                     const double tol,
                                     const double qr_tol,
                                     const double nu_max);

// logistic regression
// return the coefficients and SEs
Rcpp::List calc_coefSE_binreg(const Rcpp::NumericMatrix& X,
                              const Rcpp::NumericVector& y,
                              const int maxit,
                              const double tol,
                              const double qr_tol,
                              const double nu_max);

// logistic regression
// return (llik, individual contributions to llik, fitted probabilities, coef, SE
Rcpp::List fit_binreg(const Rcpp::NumericMatrix& X,
                      const Rcpp::NumericVector& y,
                      const bool se, // whether to include SEs
                      const int maxit,
                      const double tol,
                      const double qr_tol,
                      const double nu_max);


#endif // BINREG_H
