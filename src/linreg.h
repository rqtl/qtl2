// generic linear regression functions
#ifndef LINREG_H
#define LINREG_H

#include <RcppEigen.h>

// Calculate vector of residual sum of squares (RSS) from linear regression of Y vs X
Rcpp::NumericVector calc_rss_linreg(const Rcpp::NumericMatrix& X,
                                    const Rcpp::NumericMatrix& Y,
                                    const double tol);

// Calculate just the coefficients from linear regression of y on X
Rcpp::NumericVector calc_coef_linreg(const Rcpp::NumericMatrix& X,
                                     const Rcpp::NumericVector& y,
                                     const double tol);

// Calculate the coefficients and SEs from linear regression of y on X
Rcpp::List calc_coefSE_linreg(const Rcpp::NumericMatrix& X,
                              const Rcpp::NumericVector& y,
                              const double tol);

// Calculate matrix of residuals from linear regression of Y on X
Rcpp::NumericMatrix calc_resid_linreg(const Rcpp::NumericMatrix& X,
                                      const Rcpp::NumericMatrix& Y,
                                      const double tol);

// use calc_resid_linreg for a 3-dim array
Rcpp::NumericVector calc_resid_linreg_3d(const Rcpp::NumericMatrix& X,
                                         const Rcpp::NumericVector& P,
                                         const double tol);

// Linear regression of y on X, returning a bunch of stuff
// output is list of (coef, fitted, resid, rss, sigma, rank, df, SE)
Rcpp::List fit_linreg(const Rcpp::NumericMatrix& X,
                      const Rcpp::NumericVector& y,
                      const double tol);

#endif // LINREG_H
