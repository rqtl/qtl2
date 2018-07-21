// logistic regression via RcppEigen
#ifndef BINREG_EIGEN_H
#define BINREG_EIGEN_H

#include <RcppEigen.h>

// logistic regression by "LLt" Cholesky decomposition
// return just the log likelihood
//
// maxit = maximum number of iterations
// tol = tolerance for convergence
// eta_max = maximum value for linear predictor
double calc_ll_binreg_eigenchol(const Rcpp::NumericMatrix& X,
                                const Rcpp::NumericVector& y,
                                const int maxit,
                                const double tol,
                                const double eta_max);

// logistic regression by Qr decomposition with column pivoting
// return just the log likelihood
//
// maxit = maximum number of iterations
// tol = tolerance for convergence
// qr_tol = tolerance for QR decomposition
// eta_max = maximum value for linear predictor
//
double calc_ll_binreg_eigenqr(const Rcpp::NumericMatrix& X,
                              const Rcpp::NumericVector& y,
                              const int maxit,
                              const double tol,
                              const double qr_tol,
                              const double eta_max);

// logistic regression
// return just the coefficients
Rcpp::NumericVector calc_coef_binreg_eigenqr(const Rcpp::NumericMatrix& X,
                                             const Rcpp::NumericVector& y,
                                             const int maxit,      // max iterations
                                             const double tol,     // tolerance for convergence
                                             const double qr_tol,  // tolerance for QR decomp
                                             const double eta_max); // max value for linear predictor

// logistic regression
// return the coefficients and SEs
Rcpp::List calc_coefSE_binreg_eigenqr(const Rcpp::NumericMatrix& X,
                                      const Rcpp::NumericVector& y,
                                      const int maxit,      // max iterations
                                      const double tol,     // tolerance for convergence
                                      const double qr_tol,  // tolerance for QR decomp
                                      const double eta_max); // max value for linear predictor

// logistic regression
// return (llik, individual contributions to llik, fitted probabilities, coef, SE
Rcpp::List fit_binreg_eigenqr(const Rcpp::NumericMatrix& X,
                              const Rcpp::NumericVector& y,
                              const bool se,        // whether to include SEs
                              const int maxit,      // max iterations
                              const double tol,     // tolerance for convergence
                              const double qr_tol,  // tolerance for QR decomp
                              const double eta_max); // max value for linear predictor

#endif // BINREG_EIGEN_H
