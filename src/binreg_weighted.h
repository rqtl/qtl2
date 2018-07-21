// (weighted) logistic regression
#ifndef BINREG_WEIGHTED_H
#define BINREG_WEIGHTED_H

#include <RcppEigen.h>

// logistic regression
// return just the log likelihood
// this version uses weights
double calc_ll_binreg_weighted(const Rcpp::NumericMatrix& X,
                               const Rcpp::NumericVector& y,
                               const Rcpp::NumericVector &weights,
                               const int maxit,
                               const double tol,
                               const double qr_tol,
                               const double nu_max);

// logistic regression
// return just the coefficients
// this version uses weights
Rcpp::NumericVector calc_coef_binreg_weighted(const Rcpp::NumericMatrix& X,
                                              const Rcpp::NumericVector& y,
                                              const Rcpp::NumericVector &weights,
                                              const int maxit,
                                              const double tol,
                                              const double qr_tol,
                                              const double nu_max);

// logistic regression
// return the coefficients and SEs
// this version uses weights
Rcpp::List calc_coefSE_binreg_weighted(const Rcpp::NumericMatrix& X,
                                       const Rcpp::NumericVector& y,
                                       const Rcpp::NumericVector &weights,
                                       const int maxit,
                                       const double tol,
                                       const double qr_tol,
                                       const double nu_max);

// logistic regression
// return (llik, individual contributions to llik, fitted probabilities, coef, SE
// this version uses weights
Rcpp::List fit_binreg_weighted(const Rcpp::NumericMatrix& X,
                               const Rcpp::NumericVector& y,
                               const Rcpp::NumericVector &weights,
                               const bool se, // whether to include SEs
                               const int maxit,
                               const double tol,
                               const double qr_tol,
                               const double nu_max);

#endif // BINREG_WEIGHTED_H
