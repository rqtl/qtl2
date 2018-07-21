// logistic regression via RcppEigen
#ifndef BINREG_WEIGHTED_EIGEN_H
#define BINREG_WEIGHTED_EIGEN_H

#include <RcppEigen.h>

// logistic regression by "LLt" Cholesky decomposition
// return just the log likelihood
// this version with weights
double calc_ll_binreg_weighted_eigenchol(const Rcpp::NumericMatrix& X,
                                         const Rcpp::NumericVector& y,
                                         const Rcpp::NumericVector& weights,
                                         const int maxit,
                                         const double tol,
                                         const double eta_max);

// logistic regression by Qr decomposition with column pivoting
// return just the log likelihood
// this version with weights
double calc_ll_binreg_weighted_eigenqr(const Rcpp::NumericMatrix& X,
                                       const Rcpp::NumericVector& y,
                                       const Rcpp::NumericVector& weights,
                                       const int maxit,
                                       const double tol,
                                       const double qr_tol,
                                       const double eta_max);

// logistic regression
// return just the coefficients
// this version with weights
Rcpp::NumericVector calc_coef_binreg_weighted_eigenqr(const Rcpp::NumericMatrix& X,
                                                      const Rcpp::NumericVector& y,
                                                      const Rcpp::NumericVector& weights,
                                                      const int maxit,
                                                      const double tol,
                                                      const double qr_tol,
                                                      const double eta_max);

// logistic regression
// return the coefficients and SEs
// this version with weights
Rcpp::List calc_coefSE_binreg_weighted_eigenqr(const Rcpp::NumericMatrix& X,
                                               const Rcpp::NumericVector& y,
                                               const Rcpp::NumericVector& weights,
                                               const int maxit,
                                               const double tol,
                                               const double qr_tol,
                                               const double eta_max);

// logistic regression
// return (llik, individual contributions to llik, fitted probabilities, coef, SE
// this version with weights
Rcpp::List fit_binreg_weighted_eigenqr(const Rcpp::NumericMatrix& X,
                                       const Rcpp::NumericVector& y,
                                       const Rcpp::NumericVector& weights,
                                       const bool se, // whether to include SEs
                                       const int maxit,
                                       const double tol,
                                       const double qr_tol,
                                       const double eta_max);

#endif // BINREG_WEIGHTED_EIGEN_H
